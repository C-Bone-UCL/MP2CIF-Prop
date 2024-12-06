MP2CIF-Prop
==============

## Overview
This repository contains the scripts necessary to query the [Materials Project](https://next-gen.materialsproject.org/). With the scripts you can create a few types of datasets: CIF files with bandgap information appended to the end. .pkl tables with rows where you can find database info:

| 'Database' | 'Reduced Formula' | 'CIF' | 'Bandagap (eV)' |

The queried properties can be easily changed with a fewe script adjustements.

## Installation

To set up the environment, install the required dependencies listed in requirements.txt

```shell
pip install -r requirements.txt
```

## Usage

### Generating CIFs table

Of the format | 'Database' | 'Reduced Formula' | 'CIF' |

```shell
python gen_cif_table.py --api_key YOUR_API_KEY --size 100 --output_dir ./example_output
```

### Appending Property information to table

To get | 'Database' | 'Reduced Formula' | 'CIF' | 'Bandagap (eV)' |

```shell
python gen_cif_prop_table.py --api_key YOUR_API_KEY --size 100 --input_file ./example_table_output/materials_data.pkl --output_file ./example_table_output/CIFs_db_MP.pkl
```

### Making the CIFs extended with bandgap infomation dataset

To get: CIF + \n + Bandgap_eV 0.10 + \n + \n

```shell
python gen_cifextd.py --api_key YOUR_API_KEY --size 100 --output_dir ./example_output
```

### Dataset Generation and Analysis

The [Generate_Datasets.ipynb](https://github.com/C-Bone-UCL/MP2CIF-Prop/blob/main/Generate_Datasets.ipynb) provides a comprehensive workflow for generating and analyzing datasets. Open the notebook in Jupyter and follow the steps provided.

## Queries

For any queries feel free to email at cyprien.bone.24@ucl.ac.uk or open a GitHub issue
