
from optimade.adapters.structures.utils import (
    cell_to_cellpar,
    fractional_coordinates,
    valid_lattice_vector,
)
from optimade.models import Species as OptimadeStructureSpecies
from optimade.models import StructureResource as OptimadeStructure
import os
import pickle
import argparse
import pandas as pd
import requests
from tqdm import tqdm
from pymatgen.ext.optimade import OptimadeRester
from optimade.models import StructureResource

try:
    import numpy as np
except ImportError:
    from warnings import warn

    from optimade.adapters.warnings import AdapterPackageNotFound

    np = None  # type: ignore[assignment]
    NUMPY_NOT_FOUND = "NumPy not found, cannot convert structure to CIF"


__all__ = ("get_cif",)

# API Endpoints dictionary
api_endpoints = {
    'aflow': 'https://aflow.org/API/optimade/',
    'alexandria': 'https://alexandria.icams.rub.de/pbe/',
    'alexandria.pbe': 'https://alexandria.icams.rub.de/pbe/',
    'alexandria.pbesol': 'https://alexandria.icams.rub.de/pbesol/',
    'cod': 'https://www.crystallography.net/cod/optimade/',
    'cmr': 'https://cmr-optimade.fysik.dtu.dk/',
    'mcloud.mc3d': 'https://aiida.materialscloud.org/mc3d/optimade/',
    'mcloud.mc2d': 'https://aiida.materialscloud.org/mc2d/optimade/',
    'mcloud.2dtopo': 'https://aiida.materialscloud.org/2dtopo/optimade/',
    'mcloud.tc-applicability': 'https://aiida.materialscloud.org/tc-applicability/optimade/',
    'mcloud.pyrene-mofs': 'https://aiida.materialscloud.org/pyrene-mofs/optimade/',
    'mcloud.curated-cofs': 'https://aiida.materialscloud.org/curated-cofs/optimade/',
    'mcloud.stoceriaitf': 'https://aiida.materialscloud.org/stoceriaitf/optimade/',
    'mcloud.scdm': 'https://aiida.materialscloud.org/autowannier/optimade/',
    'mcloud.tin-antimony-sulfoiodide': 'https://aiida.materialscloud.org/tin-antimony-sulfoiodide/optimade/',
    'mcloud.optimade-sample': 'https://aiida.materialscloud.org/optimade-sample/optimade/',
    'mp': 'https://optimade.materialsproject.org/',
    'mpdd': 'http://mpddoptimade.phaseslab.org/',
    'mpds': 'https://api.mpds.io/',
    'nmd': 'https://nomad-lab.eu/prod/rae/optimade/',
    'odbx': 'https://optimade.odbx.science/',
    'odbx.odbx_misc': 'https://optimade-misc.odbx.science/',
    'odbx.gnome': 'https://optimade-gnome.odbx.science/',
    'omdb.omdb_production': 'https://optimade.openmaterialsdb.se/',
    'oqmd': 'https://oqmd.org/optimade/',
    'jarvis': 'https://jarvis.nist.gov/optimade/jarvisdft/',
    'tcod': 'https://www.crystallography.net/tcod/optimade/',
    'twodmatpedia': 'http://optimade.2dmatpedia.org/'
}

def fetch_structures(api_name, base_url, page_limit=1):
    structures = []
    endpoint = f"{base_url}/structures"
    params = {'page_limit': page_limit}
    try:
        response = requests.get(endpoint, params=params)
        response.raise_for_status()
        data = response.json()
        structures = data.get("data", [])
    except requests.exceptions.RequestException as e:
        print(f"Error fetching structures from {api_name}: {e}")
    return structures


def get_cif(
    optimade_structure: OptimadeStructure,
) -> str:
    """Get CIF file as string from OPTIMADE structure.

    Parameters:
        optimade_structure: OPTIMADE structure.

    Returns:
        The CIF file as a single Python `str` object.

    """
    # NumPy is needed for calculations
    if globals().get("np", None) is None:
        warn(NUMPY_NOT_FOUND, AdapterPackageNotFound)
        return None  # type: ignore[return-value]

    cif = """#
# Created from an OPTIMADE structure.
#
"""

    cif += f"data_{optimade_structure.id}\n\n"

    attributes = optimade_structure.attributes

    # Do this only if there's three non-zero lattice vectors
    # NOTE: This also negates handling of lattice_vectors with null/None values
    if valid_lattice_vector(attributes.lattice_vectors):  # type:ignore[arg-type]
        a_vector, b_vector, c_vector, alpha, beta, gamma = cell_to_cellpar(
            attributes.lattice_vectors  # type: ignore[arg-type]
        )

        cif += (
            f"_cell_length_a                    {a_vector:g}\n"
            f"_cell_length_b                    {b_vector:g}\n"
            f"_cell_length_c                    {c_vector:g}\n"
            f"_cell_angle_alpha                 {alpha:g}\n"
            f"_cell_angle_beta                  {beta:g}\n"
            f"_cell_angle_gamma                 {gamma:g}\n\n"
        )
        cif += (
            "_symmetry_space_group_name_H-M    'P 1'\n"
            "_symmetry_int_tables_number       1\n\n"
            "loop_\n"
            "  _symmetry_equiv_pos_as_xyz\n"
            "  'x, y, z'\n\n"
        )

        # Since some structure viewers are having issues with cartesian coordinates,
        # we calculate the fractional coordinates if this is a 3D structure and we have all the necessary information.
        if not hasattr(attributes, "fractional_site_positions"):
            attributes.fractional_site_positions = fractional_coordinates(
                cell=attributes.lattice_vectors,  # type:ignore[arg-type]
                cartesian_positions=attributes.cartesian_site_positions,  # type:ignore[arg-type]
            )

    # NOTE: This is otherwise a bit ahead of its time, since this OPTIMADE property is part of an open PR.
    # See https://github.com/Materials-Consortia/OPTIMADE/pull/206
    coord_type = (
        "fract" if hasattr(attributes, "fractional_site_positions") else "Cartn"
    )

    cif += (
        "loop_\n"
        "  _atom_site_type_symbol\n"  # species.chemical_symbols
        "  _atom_site_label\n"  # species.name + unique int
        "  _atom_site_occupancy\n"  # species.concentration
        f"  _atom_site_{coord_type}_x\n"  # cartesian_site_positions
        f"  _atom_site_{coord_type}_y\n"  # cartesian_site_positions
        f"  _atom_site_{coord_type}_z\n"  # cartesian_site_positions
        "  _atom_site_thermal_displace_type\n"  # Set to 'Biso'
        "  _atom_site_B_iso_or_equiv\n"  # Set to 1.0:f
    )

    if coord_type == "fract":
        sites = attributes.fractional_site_positions
    else:
        sites = attributes.cartesian_site_positions

    species: dict[str, OptimadeStructureSpecies] = {
        species.name: species
        for species in attributes.species  # type: ignore[union-attr]
    }

    symbol_occurences: dict[str, int] = {}
    for site_number in range(attributes.nsites):  # type: ignore[arg-type]
        species_name = attributes.species_at_sites[site_number]  # type: ignore[index]
        site = sites[site_number]

        current_species = species[species_name]

        for index, symbol in enumerate(current_species.chemical_symbols):
            if symbol == "vacancy":
                continue

            if symbol in symbol_occurences:
                symbol_occurences[symbol] += 1
            else:
                symbol_occurences[symbol] = 1
            label = f"{symbol}{symbol_occurences[symbol]}"

            cif += (
                f"  {symbol} {label} {current_species.concentration[index]:6.4f} {site[0]:8.5f}  "
                f"{site[1]:8.5f}  {site[2]:8.5f}  {'Biso':4}  {'1.000':6}\n"
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
        print(f"Error saving CIF for structure {structure_id}: {e}")

def compile_metadata(metadata_table, reduced_formula, database_name, cif_text):
    metadata_table.append({
        'reduced_formula': reduced_formula,
        'database': database_name,
        'cif_text': cif_text
    })

def save_pickle(metadata_table, filename):
    df = pd.DataFrame(metadata_table)
    if os.path.exists(filename):
        # Load existing data and append
        with open(filename, 'rb') as f:
            existing_df = pickle.load(f)
        df = pd.concat([existing_df, df], ignore_index=True)
    try:
        with open(filename, 'wb') as f:
            pickle.dump(df, f)
        print(f"Metadata successfully saved to {filename}")
    except Exception as e:
        print(f"Error saving pickle file: {e}")

def main():
    parser = argparse.ArgumentParser(description="Fetch structures from OPTIMADE databases and convert to CIF.")
    parser.add_argument('--database', nargs='+', required=True, help="List of databases to fetch data from.")
    parser.add_argument('--dataset', type=int, required=True, help="Number of entries to fetch from each database.")
    parser.add_argument('--output_dir', required=True, help="Directory to save CIF files and metadata pickle file.")
    args = parser.parse_args()

    output_cif_dir = os.path.join(args.output_dir, 'cif_files')
    output_pkl_file = os.path.join(args.output_dir, 'structures_metadata.pkl')

    metadata_table = []

    for db in args.database:
        if db not in api_endpoints:
            print(f"Database '{db}' not recognized. Skipping.")
            continue

        base_url = api_endpoints[db]
        print(f"Fetching structures from {db}...")
        structures = fetch_structures(db, base_url, page_limit=args.dataset)
        print(f"Number of structures fetched from {db}: {len(structures)}")

        for structure_dict in tqdm(structures, desc=f"Processing {db}", leave=False):
            try:
                structure = StructureResource(**structure_dict)
            except Exception as e:
                print(f"Error parsing structure: {e}")
                continue

            structure_id = structure.id
            attributes = structure.attributes
            reduced_formula = attributes.chemical_formula_reduced or 'N/A'

            try:
                cif_text = get_cif(structure)
            except Exception as e:
                print(f"Error generating CIF for structure {structure_id}: {e}")
                cif_text = None

            if cif_text:
                try:
                    save_cif(cif_text, output_cif_dir, structure_id)
                    compile_metadata(metadata_table, reduced_formula, db, cif_text)
                except Exception as e:
                    print(f"Error saving CIF or compiling metadata for structure {structure_id}: {e}")
            else:
                print(f"Skipping structure {structure_id} due to conversion issues.")

    # Save all metadata to a pickle file
    try:
        save_pickle(metadata_table, output_pkl_file)
    except Exception as e:
        print(f"Error saving metadata: {e}")

if __name__ == "__main__":
    main()
