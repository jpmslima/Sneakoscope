# Sneakoscope: README

The **Sneakoscope** is a comprehensive Streamlit-based web application designed for protein variant analysis, pathogenicity prediction, and 3D structural modeling. It integrates data from UniProt, AlphaMissense, and SWISS-MODEL to provide an interactive dashboard for researchers to explore how specific mutations might affect protein structure and function.

---

## Key Features

### 1. Data Integration & Visualization
* **Multi-Source Fetching:** Automatically retrieves protein features and variation data from the EBI Proteins API.
* **Mutation Frequency:** Generates a bar chart showing the frequency of known mutations at every position along the protein sequence.
* **Interactive Feature Map:** A staggered, multi-lane visualization of protein domains, sites, and regions (e.g., active sites, binding domains).
* **AlphaMissense Pathogenicity:** Plots predicted pathogenicity scores for variants, categorizing them as "likely benign," "ambiguous," or "likely pathogenic".

### 2. Mutation Constructor & Alignment
* **Sequence Engineering:** Allows users to select multiple known mutations and apply them to the wild-type reference sequence.
* **Visual Alignment:** Displays a side-by-side comparison of the original and mutated sequences, highlighting changes in red for easy identification.
* **FASTA Export:** Generates a ready-to-use FASTA format string for the mutated protein.

### 3. 3D Structure Modeling & Quality Control
* **SWISS-MODEL Integration:** Submits the mutated sequence to the SWISS-MODEL API to generate a new 3D structure based on AlphaFold templates.
* **Interactive 3D Comparison:** Uses `st-molstar` to provide side-by-side synchronized views of the wild-type and mutated structures, colored by pLDDT confidence scores.
* **Ramachandran Analysis:** Computes and plots $\phi$ (Phi) and $\psi$ (Psi) torsion angles to assess the stereochemical quality of the generated models.


---

## Technical Architecture

### Core Libraries
| Library | Purpose |
| :--- | :--- |
| **Streamlit** | Web interface and state management. |
| **Pandas / NumPy** | Data manipulation and numerical calculations. |
| **Plotly** | Interactive graphing (Bar, Scatter, and Feature Maps). |
| **Bio.PDB** | Parsing PDB files and calculating torsion angles. |
| **st-molstar** | High-performance 3D molecular visualization. |

### Data Flow
1.  **Input:** User provides a UniProt Accession ID (e.g., `P63000`).
2.  **Collection:** The app fetches JSON features, variation data, and AlphaMissense CSVs.
3.  **Modeling:** If the user provides a SWISS-MODEL API token, the app initiates a remote modeling job using the mutated sequence.
4.  **Analysis:** The app calculates B-factors (pLDDT) and Ramachandran plots for both the wild-type and the newly generated mutant.

---

## Getting Started

### Prerequisites
* Python 3.8+
* A SWISS-MODEL API Token (optional, for 3D modeling).

### Installation
```bash
pip install streamlit pandas requests plotly biopython streamlit-molstar
```

### Running the App
1. Save the code as `app.py`.
2. Run the command: `streamlit run app.py`.
3. Enter a **UniProt ID** and click **"Cast Revelio!"** to begin.

---

## Important Notes
* **Caching:** The app uses `@st.cache_data` to ensure that repeated queries for the same protein ID are instant and do not overload external APIs.
* **pLDDT Coloring:** Structure visualizations use a standard color scale where Blue indicates high confidence ($>90$) and Red indicates low confidence ($<50$).
* **Asynchronous Jobs:** SWISS-MODEL jobs are polled every 10 seconds until completion or failure.

---

## ⚖️ License
* This project is licensed under the **MIT License**. You are free to use, modify, and distribute this software for academic or commercial purposes, provided that proper credit is given to the original author.
