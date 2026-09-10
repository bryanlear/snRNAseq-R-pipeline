## Morabito et al., 2021

Neurodegenerative disorders such as **AD** are characterized by **massive neuronal loss** with *gliosis*. 

* **Glial Cells**: (Neuroglia) are non-neuronal cells in the CNS and PNS that support, protect, and regulate neurons:
  * Astrocytes (CNS)
  * Microglia (CNS)
  * Oligodendrocytes (CNS)
  * Ependymal Cells (CNS)
  * Schwann cells (PNS)
  * Satellite cells (PNS)

* **Gliosis**: involves proliferation, hypertrophy, morphological transformation of mostly astrocytes, microglia, and oligodendrocytes in response to CNS injury, infection, ischemia, neurodegenerative disease, or trauma.

GWAS studies have shown variants in *distal regulatory elements* - e.g., enhancers (often cell-type specific regions in relevant tissues) and not in genes. 

**GWAS** identifies DNA variants associated with a trait or disease. In order to determine the cell type affected by $x$ variant, one must know:

- Regulatory element containing the variant 
- Gene that controls the element
- Which cell types the regulatory *connection* is active

Standard GWAS cannot show the cell type that uses the regulatory element or the gene the element regulates.

Therefore, one needs a map that connects distant regulatory elements to target genes. ([FUMA](https://fuma.ctglab.nl/tutorial)) can be suitable for an initial GWAS to gene and cell type analysis (it uses positional + eQTL + chromatin interaction data) but at higher resolution (e.g., cell type- specific enhancer-gene map), snATAC-seq can identify accessible regulatory elements in each cell type.

**snATAC-seq** shows where chromatin is accessible in each cell type which helps identify cell type in which a GWAS variant can have a regulatory activity. 

Thus, **co-accessibility** analysis and integration with **snRNA-seq** can provide stronger evidence for candidate ***target*** gene.

* **Co-accessibility**: test whether 2 DNA regions tend to be open at the same time across many cells. (e.g., A= distant enhancer, B= gene promoter) --> Correlated accessibility would imply enhancer and promoter may be part of same regulatory system. 

---

*Isolated nuclei from prefontal cortex*: 

- snATAC-seq --> late-stage AD n=12 and control n=7
- snRNA-seq --> late-stage AD n=11 and control n=7

**Braak staging**: measures anatomical spread to *tau neurofibrillary tangles**:
  - Stages I-II: mainly transentorhinal regions
  - Stages III-IV: spread into limbic regions
  - Stages V-VI: extensuve neocortical involvement
  
**Plaque staging**: measures amount/anatomical spread of amyloid-beta plaques. Higher stage/score = greater plaque burden (systems include Thal phases and CERAD plaque scores).

**Gene program/metagene**: group of genes that changes together

---

snATAC-seq matrix --> *nuclei x accessible DNA regions*:

* `1`: region was detected as accessible
* `0`: no accessibility was detected

**LIS**:
1. TF-IDF (comes from text analysis)
   - A region accessible in almost every nucleus gets less weight
   - A region accessible in a limited set of nuclei receives more weight
   - Also differences in total number of fragments per nucleus are reduced
2. SVD
   - Compression of weighted region matrix into small number of LSI components:
      
    130k nuclei x 200k regions Matrix $\rightarrow$ 130k x 30-50 LSI components

    each component = pattern chromatin accessibility 
3. Morabito et al., applies **MNN** to LSI matrix: MNN searches for cells in different batches that are each other's nearest neighbors. Pairs likely represent same cell type. MNN the uses differences between matched cells to estimate local batch correction vectors. 

---

