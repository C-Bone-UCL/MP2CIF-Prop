import argparse
import os
import time
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter
from requests.exceptions import HTTPError
from tqdm import tqdm

def generate_and_save_data(api_key, num_samples, output_dir):
    mpr = MPRester(api_key)
    materials = mpr.query({}, ['material_id', 'pretty_formula'])

    if num_samples != 'max':
        materials = materials[:int(num_samples)]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    

    for material in tqdm(materials):
        material_id = material["material_id"]
        try:
            structure = mpr.get_structure_by_material_id(material_id)
            bandgap_info = mpr.query(criteria={"material_id": material_id}, properties=["band_gap", "pretty_formula"])
            time.sleep(0.1)  # Add a delay to avoid API rate limiting
        except Exception as e:  # Catch HTTPError and other possible exceptions
            print(f"Error retrieving data for material {material_id}: {e}")
            continue  # Skip this material

        bandgap = bandgap_info[0].get("band_gap", "N/A")
        bandgap = round(bandgap, 3) if isinstance(bandgap, (float, int)) else bandgap

        filename = os.path.join(output_dir, f"{material['pretty_formula']}_{material_id}.cif")
        with open(filename, "w") as f:
            cif_writer = CifWriter(structure)
            cif_data = cif_writer.__str__()

            header = f"# MP Entry {material_id} {material['pretty_formula']}\n"
            cif_data = header + cif_data + f"\nBandgap_eV {bandgap}\n"

            f.write(cif_data)

def main():
    parser = argparse.ArgumentParser(description="Generate dataset from Materials Project and save as CIF files.")
    parser.add_argument("--api_key", type=str, required=True, 
                        help="API key for accessing the Materials Project database.")
    parser.add_argument("--size", type=str, default="100", 
                        help="Specify the dataset size or 'max' for all data. Default is 100.")
    parser.add_argument("--output_dir", type=str, default="/home/uccacbo/CrystaLLM/BG_CIFs", 
                        help="Directory to save the CIF files.")
    args = parser.parse_args()

    # Determine dataset size
    if args.size.lower() == "max":
        num_samples = 'max'
    else:
        try:
            num_samples = int(args.size)
        except ValueError:
            raise ValueError("Invalid value for size. Please specify an integer or 'max'.")

    # Generate and save dataset
    generate_and_save_data(args.api_key, num_samples, args.output_dir)

if __name__ == "__main__":
    main()
