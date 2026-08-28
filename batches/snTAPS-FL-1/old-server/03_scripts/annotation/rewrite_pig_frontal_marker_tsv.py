from pathlib import Path
import pandas as pd

out = Path("03_scripts/annotation/pig_frontal_lobe_marker_symbols.tsv")

markers = [
    ("Excitatory_neuron","RBFOX3"),
    ("Excitatory_neuron","SNAP25"),
    ("Excitatory_neuron","SYT1"),
    ("Excitatory_neuron","MAP2"),
    ("Excitatory_neuron","SLC17A7"),
    ("Excitatory_neuron","CAMK2A"),
    ("Excitatory_neuron","GRIN1"),
    ("Excitatory_neuron","NEFL"),
    ("Excitatory_neuron","NEFM"),

    ("Upper_layer_excitatory","SATB2"),
    ("Upper_layer_excitatory","CUX1"),
    ("Upper_layer_excitatory","CUX2"),
    ("Upper_layer_excitatory","RORB"),

    ("Deep_layer_excitatory","BCL11B"),
    ("Deep_layer_excitatory","FEZF2"),
    ("Deep_layer_excitatory","TBR1"),
    ("Deep_layer_excitatory","FOXP2"),

    ("Inhibitory_neuron","GAD1"),
    ("Inhibitory_neuron","GAD2"),
    ("Inhibitory_neuron","SLC32A1"),
    ("Inhibitory_neuron","PVALB"),
    ("Inhibitory_neuron","SST"),
    ("Inhibitory_neuron","VIP"),
    ("Inhibitory_neuron","LAMP5"),
    ("Inhibitory_neuron","RELN"),
    ("Inhibitory_neuron","NPY"),

    ("Astrocyte","AQP4"),
    ("Astrocyte","GFAP"),
    ("Astrocyte","ALDH1L1"),
    ("Astrocyte","SLC1A2"),
    ("Astrocyte","SLC1A3"),
    ("Astrocyte","GJA1"),
    ("Astrocyte","SOX9"),

    ("Oligodendrocyte","MBP"),
    ("Oligodendrocyte","MOG"),
    ("Oligodendrocyte","PLP1"),
    ("Oligodendrocyte","MAG"),
    ("Oligodendrocyte","MOBP"),
    ("Oligodendrocyte","CLDN11"),

    ("OPC","PDGFRA"),
    ("OPC","CSPG4"),
    ("OPC","VCAN"),
    ("OPC","OLIG1"),
    ("OPC","OLIG2"),
    ("OPC","SOX10"),

    ("Microglia","AIF1"),
    ("Microglia","CX3CR1"),
    ("Microglia","P2RY12"),
    ("Microglia","CSF1R"),
    ("Microglia","TYROBP"),
    ("Microglia","TREM2"),
    ("Microglia","C1QA"),
    ("Microglia","C1QB"),
    ("Microglia","C1QC"),

    ("Endothelial","PECAM1"),
    ("Endothelial","CLDN5"),
    ("Endothelial","VWF"),
    ("Endothelial","KDR"),
    ("Endothelial","FLT1"),

    ("Pericyte_SMC","PDGFRB"),
    ("Pericyte_SMC","RGS5"),
    ("Pericyte_SMC","ACTA2"),
    ("Pericyte_SMC","TAGLN"),
]

df = pd.DataFrame(markers, columns=["group", "gene"])
df.to_csv(out, sep="\t", index=False)
print("[INFO] wrote", out)
print(df.head())
