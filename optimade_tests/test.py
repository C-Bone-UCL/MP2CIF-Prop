from optimade.client import OptimadeClient
from optimade.models import Species as OptimadeStructureSpecies
from optimade.models import StructureResource as OptimadeStructure
from pymatgen.ext.optimade import OptimadeRester
from optimade.adapters.structures.utils import (
    cell_to_cellpar,
    fractional_coordinates,
    valid_lattice_vector,
)
import os
import pickle
import argparse
import pandas as pd
from tqdm import tqdm
from warnings import warn
import logging

try:
    import numpy as np
except ImportError:
    from optimade.adapters.warnings import AdapterPackageNotFound

    np = None
    NUMPY_NOT_FOUND = "NumPy not found, cannot convert structure to CIF"


def log_error(url, results):
    if results.get("errors"):
        logging.error(f"Errors encountered for {url}: {results['errors']}")


def get_cif(
    optimade_structure: OptimadeStructure,
) -> str:
    """Generate CIF file as string from OPTIMADE structure with added fields."""

    if globals().get("np", None) is None:
        warn(NUMPY_NOT_FOUND, AdapterPackageNotFound)
        return None

    cif = """#
# Created from an OPTIMADE structure.
#
"""

    cif += f"data_{optimade_structure.id}\n"

    attributes = optimade_structure.attributes

    # Check for valid lattice vectors and extract cell parameters
    if valid_lattice_vector(attributes.lattice_vectors):
        a_vector, b_vector, c_vector, alpha, beta, gamma = cell_to_cellpar(
            attributes.lattice_vectors
        )
        cif += (
            f"_cell_length_a {a_vector:g}\n"
            f"_cell_length_b {b_vector:g}\n"
            f"_cell_length_c {c_vector:g}\n"
            f"_cell_angle_alpha {alpha:g}\n"
            f"_cell_angle_beta {beta:g}\n"
            f"_cell_angle_gamma {gamma:g}\n"
        )
    
    # Space group and symmetry fields
    if attributes.space_group_symbol_hermann_mauguin:
        cif += f"_symmetry_space_group_name_H-M '{attributes.space_group_symbol_hermann_mauguin}'\n"
    if attributes.space_group_it_number:
        cif += f"_symmetry_int_tables_number {attributes.space_group_it_number}\n"

    # Basic chemical formula and metadata
    cif += (
        f"_chemical_formula_structural {attributes.chemical_formula_descriptive}\n"
        f"_chemical_formula_sum '{attributes.chemical_formula_reduced}'\n"
    )

    # Cell volume and units (Z)
    if hasattr(attributes, "cell_volume"):
        cif += f"_cell_volume {attributes.cell_volume:g}\n"
    if hasattr(attributes, "cell_formula_units_Z"):
        cif += f"_cell_formula_units_Z {attributes.cell_formula_units_Z}\n"

    # Symmetry operations (if provided)
    if attributes.space_group_symmetry_operations_xyz:
        cif += "loop_\n"
        cif += "_symmetry_equiv_pos_site_id\n"
        cif += "_symmetry_equiv_pos_as_xyz\n"
        for i, operation in enumerate(attributes.space_group_symmetry_operations_xyz, start=1):
            cif += f" {i} '{operation}'\n"

    cif += "loop_\n"
    cif += "_atom_site_type_symbol\n"
    cif += "_atom_site_label\n"
    cif += "_atom_site_symmetry_multiplicity\n"
    cif += "_atom_site_occupancy\n"

    coord_type = (
        "fract" if hasattr(attributes, "fractional_site_positions") else "Cartn"
    )
    cif += (
        f"_atom_site_{coord_type}_x\n"
        f"_atom_site_{coord_type}_y\n"
        f"_atom_site_{coord_type}_z\n"
    )

    # Extract fractional or Cartesian coordinates as appropriate
    sites = (
        attributes.fractional_site_positions
        if hasattr(attributes, "fractional_site_positions")
        else attributes.cartesian_site_positions
    )

    # Extract species data
    species_dict = {species.name: species for species in attributes.species}

    symbol_occurrences = {}
    for site_number in range(attributes.nsites):
        species_name = attributes.species_at_sites[site_number]
        site = sites[site_number]
        current_species = species_dict[species_name]

        for idx, symbol in enumerate(current_species.chemical_symbols):
            if symbol == "vacancy":
                continue

            # Track occurrences for unique labels
            symbol_occurrences[symbol] = symbol_occurrences.get(symbol, 0) + 1
            label = f"{symbol}{symbol_occurrences[symbol]}"

            cif += (
                f"  {symbol} {label} {current_species.concentration[idx]:6.4f} "
                f"{site[0]:8.5f}  {site[1]:8.5f}  {site[2]:8.5f}\n"
            )

    return cif


def save_cif(cif_text, output_dir, structure_id):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_path = os.path.join(output_dir, f"{structure_id}.cif")
    try:
        with open(file_path, 'w') as f:
            f.write(cif_text)
    except Exception as e:
        logging.error(f"Error saving CIF for structure {structure_id}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch structures from OPTIMADE databases and convert to CIF."
    )
    parser.add_argument(
        '--dataset',
        type=int,
        required=True,
        help="Number of entries to fetch from each database.",
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help="Directory to save CIF files and metadata pickle file.",
    )
    args = parser.parse_args()

    # make sure directory exists
    # os.mkdir -p args.output_dir

    output_cif_dir = os.path.join(args.output_dir, 'cif_files')
    output_pkl_file = os.path.join(args.output_dir, 'structures_metadata.pkl')

    log_file = os.path.join(args.output_dir, 'outputs.log')
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    metadata_table = []

    # Essential response fields for additional details
    response_fields = [
        'id', 'lattice_vectors', 'cartesian_site_positions', 'species',
        'species_at_sites', 'nsites', 'last_modified', 'structure_features',
        'elements', 'chemical_formula_descriptive', 'chemical_formula_reduced',
        'cell_volume', 'cell_formula_units_Z', 'space_group_symbol_hermann_mauguin',
        'space_group_it_number', 'space_group_symmetry_operations_xyz',
    ]

    client = OptimadeClient(
        max_results_per_provider=args.dataset,
        callbacks=[log_error],
        use_async=True,
        silent=False,
        max_attempts=3,
        exclude_providers=['tcod','cod']
    )

    try:
        all_results = client.get(response_fields=response_fields)
        logging.info("Successfully fetched structures from OPTIMADE providers.")
    except Exception as e:
        logging.error(f"Error fetching structures: {e}")
        # handle the ones that make an error
        # try:
        #     # log which responses_fields failed, then requery with only the ones that were successful
        #     response_fields = [field for field in response_fields if field not in all_results['errors']]
        #     all_results = client.get(response_fields=response_fields)
        #     logging.info("Successfully fetched structures from OPTIMADE providers.")
        # except Exception as e:
        #     logging.error(f"Error fetching structures: {e}")

        return

    structures_results = all_results.get('structures', {})
    total_structures = 0

    for filter_string, filter_results in structures_results.items():
        for base_url, results in tqdm(
            filter_results.items(), desc="Processing databases"
        ):
            if "errors" in results and results["errors"]:
                logging.error(f"Errors encountered for {base_url}: {results['errors']}")
                
            data = results.get("data", [])
            for structure_dict in data:
                try:
                    structure = OptimadeStructure(**structure_dict)
                except Exception as e:
                    logging.error(f"Error parsing structure: {e}")
                    continue

                structure_id = structure.id
                reduced_formula = (
                    structure.attributes.chemical_formula_reduced or 'N/A'
                )
                logging.info(f"Processing structure {structure_id} with formula {reduced_formula}")

                try:
                    cif_text = get_cif(structure)
                    if cif_text:
                        save_cif(cif_text, output_cif_dir, structure_id)
                        metadata_table.append(
                            {
                                'reduced_formula': reduced_formula,
                                'database': base_url,
                                'cif_text': cif_text,
                            }
                        )
                        logging.info(f"Successfully processed structure {structure_id} from {base_url}")
                        total_structures += 1
                except Exception as e:
                    logging.error(f"Error handling structure {structure_id}: {e}")

    if metadata_table:
        df = pd.DataFrame(metadata_table)
        try:
            with open(output_pkl_file, 'wb') as f:
                pickle.dump(df, f)
            logging.info(f"Metadata successfully saved to {output_pkl_file}")
        except Exception as e:
            logging.error(f"Error saving metadata: {e}")
    else:
        logging.warning("No metadata to save.")

    logging.info(f"Total CIF files generated: {total_structures}")

if __name__ == '__main__':
    main()
