
import streamlit as st
from Bio.SeqUtils import molecular_weight
from collections import Counter

# Page configuration
st.set_page_config(
    page_title="Protein Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

# Application title
st.title("🧬 Protein Sequence Analyzer")

st.subheader("A Bioinformatics tool for protein sequence analysis")

st.write("Welcome! Analyze your protein sequence using Python and Biopython.")

# Sequence input
st.header("Enter Your Protein Sequence")

sequence = st.text_area(
    "Paste your amino acid sequence below:",
    placeholder="Example: MKWVTFISLLFLFSSAYS",
    height=150
)

# Analyze button
if st.button("Analyze Protein Sequence"):

    if sequence.strip():

        protein = sequence.strip().upper().replace(" ", "").replace("\n", "")

        # Validate sequence
        valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

        if not all(aa in valid_amino_acids for aa in protein):

            st.error("Invalid protein sequence! Please enter valid amino acids.")

        else:

            st.success("Protein sequence analyzed successfully!")

            # 1. Display sequence
            st.subheader("Input Protein Sequence")
            st.code(protein)

            # 2. Sequence length
            length = len(protein)

            # 3. Amino acid count
            amino_acid_count = Counter(protein)

            # 4. Unique amino acids
            unique_amino_acids = len(amino_acid_count)

            # 5. Molecular weight
            weight = molecular_weight(protein, seq_type="protein")

            st.divider()

            st.subheader("Protein Analysis Results")

            col1, col2, col3 = st.columns(3)

            col1.metric("Sequence Length", length)
            col2.metric("Unique Amino Acids", unique_amino_acids)
            col3.metric("Molecular Weight", f"{weight:.2f} Da")


            # GRAVY Score Calculation

            st.divider()

            st.subheader("Hydrophobicity Analysis (GRAVY)")

            from Bio.SeqUtils.ProtParam import ProteinAnalysis

            protein_analysis = ProteinAnalysis(protein)

            gravy_score = protein_analysis.gravy()

            st.metric("GRAVY Score", f"{gravy_score:.3f}")
            
            if gravy_score > 0:
               st.success("The protein is generally hydrophobic.")

            elif gravy_score < 0:
               st.info("The protein is generally hydrophilic.")

            else:
               st.warning("The protein has a neutral average hydropathicity.")

            
            # Isoelectric Point Calculation

            st.divider()

            st.subheader("Isoelectric Point Analysis")

            isoelectric_point = protein_analysis.isoelectric_point()

            st.metric(
                "Isoelectric Point (pI)",
               f"{isoelectric_point:.2f}"
            )

            st.write(
                "The isoelectric point represents the pH "
               "at which the protein has zero net charge."
            )


            # Aromaticity Calculation

            st.divider()

            st.subheader("Aromaticity Analysis")

            aromaticity = protein_analysis.aromaticity()

            st.metric(
                "Aromaticity",
               f"{aromaticity:.4f}"
            )

            st.write(
                f"Aromatic amino acid percentage: "
                f"{aromaticity * 100:.2f}%"
            )


            # Instability Index Calculation

            st.divider()

            st.subheader("Protein Stability Analysis")

            instability_index = protein_analysis.instability_index()

            st.metric(
                "Instability Index",
               f"{instability_index:.2f}"
            )

            if instability_index < 40:
                st.success("Predicted Stable Protein")

            else:
                st.warning("Predicted Unstable Protein")

            st.write(
                "An instability index below 40 indicates a "
                "predicted stable protein, while a value above "
                "40 indicates a predicted unstable protein."
            )


            # Aliphatic Index Calculation

            st.divider()

            st.subheader("Aliphatic Index Analysis")

            A = protein.count("A") / length * 100
            V = protein.count("V") / length * 100
            I = protein.count("I") / length * 100
            L = protein.count("L") / length * 100

            aliphatic_index = A + 2.9 * V + 3.9 * (I + L)

            st.metric(
                "Aliphatic Index",
               f"{aliphatic_index:.2f}"
            )

            st.write(
                "This value estimates the relative volume "
                "of aliphatic amino acid side chains."
            )


            # Download Protein Analysis Report

            st.divider()

            st.subheader("Download Analysis Report")

            # Create report content
            report = f"""
            ========================================
                  PROTEIN SEQUENCE ANALYSIS
            ========================================

            INPUT PROTEIN SEQUENCE:
            {protein}

            ----------------------------------------
            PROTEIN PROPERTIES
            ----------------------------------------

            Sequence Length: {length} amino acids

            Molecular Weight: {weight:.2f} Da

            GRAVY Score: {gravy_score:.3f}

            Isoelectric Point (pI): {isoelectric_point:.2f}

            Aromaticity {aromaticity:.4f}

            Instability Index: {instability_index:.2f}

            Aliphatic Index: {aliphatic_index:.2f}

            ----------------------------------------
            AMINO ACID COMPOSITION
            ----------------------------------------
            """

            for amino_acid, count in amino_acid_count.items():

                percentage = (count / length) * 100

            report += (
                f"\n{amino_acid}: {count} "
                f"({percentage:.2f}%)"
            )

            report += "\n\n========================================\n"
            report += "       END OF ANALYSIS REPORT\n"
            report += "========================================\n"

            # Download button
            st.download_button(
                label="Download Protein Analysis Report",
                data=report,
                file_name="protein_analysis_report.txt",
                mime="text/plain"
            )

   
            # Amino acid composition
            st.subheader("Amino Acid Composition")

            for amino_acid, count in amino_acid_count.items():

                 percentage = (count / length) * 100

            st.write(
               f"{amino_acid} : {count} "
               f"({percentage:.2f}%)"
            )


            
            # Amino Acid Composition Chart

            st.subheader("Amino Acid Composition Chart")

            chart_data = {}

            for amino_acid, count in amino_acid_count.items():
                percentage = (count / length) * 100
                chart_data[amino_acid] = percentage

            st.write("Chart data:", chart_data)

            st.bar_chart(
               chart_data,
               x_label="Amino Acid",
               y_label="Percentage (%)"
            )



    else:

        st.warning("Please enter a protein sequence.")

st.divider()

st.caption("Developed using Python, Biopython and Streamlit")

