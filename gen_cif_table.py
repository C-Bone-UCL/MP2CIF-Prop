import argparse
import os
import time
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifWriter
from requests.exceptions import HTTPError
from tqdm import tqdm
import pandas as pd

def generate_and_save_data(api_key, num_samples, output_dir):
    mpr = MPRester(api_key)
    materials = mpr.query({}, ['material_id', 'pretty_formula'])

    if num_samples != 'max':
        materials = materials[:int(num_samples)]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    data = []  # List to store data for DataFrame

    for material in tqdm(materials):
        material_id = material["material_id"]
        try:
            structure = mpr.get_structure_by_material_id(material_id)
            time.sleep(0.1)  # Add a delay to avoid API rate limiting
        except Exception as e:
            print(f"Error retrieving data for material {material_id}: {e}")
            continue

        # Generate CIF content
        cif_writer = CifWriter(structure)
        cif_data = cif_writer.__str__()

        header = f"# MP Entry {material_id} {material['pretty_formula']}\n"
        cif_data = header + cif_data  # Exclude Bandgap info

        # Append data to list
        data.append({
            "Database": "MP",
            "Reduced Formula": material["pretty_formula"],
            "CIF": cif_data
        })

    # Save data to .pkl file
    df = pd.DataFrame(data)
    output_file = os.path.join(output_dir, "materials_data.pkl")
    df.to_pickle(output_file)
    print(f"Data saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate dataset from Materials Project and save as a .pkl file.")
    parser.add_argument("--api_key", type=str, required=True, 
                        help="API key for accessing the Materials Project database.")
    parser.add_argument("--size", type=str, default="100", 
                        help="Specify the dataset size or 'max' for all data. Default is 100.")
    parser.add_argument("--output_dir", type=str, default="/home/uccacbo/CrystaLLM/BG_CIFs", 
                        help="Directory to save the .pkl file.")
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
