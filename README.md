# `hydrogen_applications`

`hydrogen_applications` is a repository that enables the quantification of the climate-effectiveness of planned hydrogen projects and applications. It uses data obtained from the International Energy Agency (IEA) and code developed during the TRANSIENCE project (see [Acknowledgements](#acknowledgements)).

## 📁 Repository Structure

This repository contains the following files and folders:

- `0_get_renewables_battery_ratios.ipynb`  
  → Jupyter notebook to optimally design low-carbon (solar PV or wind) electrolytic hydrogen production systems.

- `1_set_up_full_environment.ipynb`  
  → Jupyter notebook to set up the (prospective) LCA databases needed during the next steps.

- `2_start_gen_db_calc_impacts.ipynb`  
  → **Main script** used to modify background LCA databases, import IEA Excel data, create life cycle inventories per facility, and generate results and figures.

- `3_1_create_additional_prospective_regionalized_dbs.ipynb`  
  → Script to generate additonal background LCA databases based on another IAM (IMAGE).

- `3_2_create_additional_figs_IMAGE.ipynb`  
  → Script to generate additonal figures based on another IAM (IMAGE).

- `3_3_create_sensitivity_figs.ipynb`  
  → Script to generate additonal figures for the sensitivity analysis considering different sector transformations in premise and for the Sankey diagrams.

- `config.py`  
  → Configuration variables.

- `functions_db.py`  
  → Functions used within the Jupyter notebooks for database generation and impact calculation.

- `functions_profile_gen.py`  
  → Functions for renewable potential calculations used in `0_get_renewables_battery_ratios.ipynb`.

- `regionalization.py`
  → Code to regionalize life cycle inventories, adapted for the purpose of this work but based on the wurst Python package.

- `mappings.py`  
  → Mappings used/shared in the different notebooks.

- `figs/`  
  → Folder containing high-quality figures from the analysis.

- `data/`  
  → Folder with input data, one needs to download the IEA hydrogen projects file (year 2024 used) via https://www.iea.org/data-and-statistics/data-product/hydrogen-production-and-infrastructure-projects-database. In this repository, the Excel file from the IEA is called:  
  `IEA Hydrogen Production Projects Database_2024.xlsx`
  → Note that the LCI files `H2-DRI_LCI.xlsx` and `BF-BOF-CCS_Carina.xlsx` were confidential at the time the paper was prepared, but that those LCI Excel files are now available (as such, those links should be replaced) via the following links based on the work of [Harpprecht et al. (2025)](https://pubs.rsc.org/en/content/articlelanding/2025/ee/d5ee01356a): 
    * [Code and data for publication: Future Environmental Impacts of Global Iron and Steel Production](https://zenodo.org/records/14968094).


Further documentation is embedded within the individual scripts. The expected outputs, generated data, and visualizations are showcased in the notebooks.

---

## 🔧 Dependencies

To run this repository or use the code, the following credentials and tools are required:

- A background [`ecoinvent database`](https://ecoinvent.org/database/), here v3.10 has been used.
- `KEY_PREMISE`: A key for [`premise`](https://github.com/polca/premise) to generate prospective LCA databases, including additional life cycle inventories.  
- `USER_PW`: An account and password to access the background LCA database **ecoinvent**.  
- A valid **Gurobi license key** to optimize low-carbon hydrogen production systems. [Setup instructions](https://support.gurobi.com/hc/en-us/articles/12872879801105-How-do-I-retrieve-and-set-up-a-Gurobi-license)  
- IEA hydrogen projects file (year 2024 used) via https://www.iea.org/data-and-statistics/data-product/hydrogen-production-and-infrastructure-projects-database.
- All required Python packages are listed in the environment file `bw_env_hydrogen.yml`.

---

## 📄 License, Citing, and Scientific References

If you use this repository, the data, or any of the included code, please cite the following paper:  
*Terlouw, T., Moretti, C., Harpprecht, C., Sacchi, R., McKenna, R., & Bauer, C. (2025). Global greenhouse gas emissions mitigation potential of existing and planned hydrogen projects. Nature Energy.*

This repository includes material derived from International Energy Agency (IEA) sources. In line with the IEA’s Creative Commons license, the authors acknowledge that they are solely responsible for this derived work, and it is not endorsed by the IEA.

For licensing information, see the `LICENSE` file.

---

## 🤝 Contributing

Contributions are welcome!  
For major suggestions, collaborations, or structural changes, please contact:

**Tom Terlouw**  
📧 [tom.terlouw@psi.ch](mailto:tom.terlouw@psi.ch)

---

## 🙏 Acknowledgements

This repository builds on several scientific contributions (see [License, citing, and scientific references](#-license-citing-and-scientific-references)) and was developed under the **TRANSIENCE** project:  
🔗 https://www.transience.eu/

It also builds upon:
- The **`premise`** framework (Sacchi et al., 2022), which enables the modification of background LCA databases.

### Supported by:

- **SHELTERED**, funded by the Swiss Federal Office of Energy (SFOE) 🇨🇭  
- **TRANSIENCE**, funded by:  
  - The European Health and Digital Executive Agency (HADEA) 🇪🇺  
  - The Swiss State Secretariat for Education, Research and Innovation (SERI) 🇨🇭  
  - The UK Research and Innovation (UKRI) Horizon Europe Guarantee 🇬🇧  
- **reFuel.ch**, funded under the SWEET programme (Grant No. SI/50271) 🇨🇭

> *The views expressed are those of the authors and do not necessarily reflect those of the European Commission or other funding institutions.*
