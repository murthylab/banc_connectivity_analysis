# **Effective Connectome Alignment For Cell Typing and Identifying Sex Differences**

Arie Matsliah1\*, Christopher K Salmon1\*, Alexander Shakeel Bates, Helen H Yang, Jasper S Phelps, Eric Perlman1, Kevin Delgado1, Benjamin Silverman1, Jay Gager1, Szi-Chieh Yu, Kyle P Willie1, Austin T Burke1, Ryan Willie1, Doug Bland1, Marissa Sorek1, Celia David1, Amy R Sterling1, Rachel I Wilson, Wei-Chung Allen Lee, H Sebastian Seung1, Mala Murthy1† and The FlyWire Consortium

1Princeton Neuroscience Institute, Princeton University, Princeton, NJ, USA  
\*These authors contributed equally  
†Corresponding authors  
e-mail: arie@princeton.edu, mmurthy@princeton.edu

Interpreting the increasing number connectomes requires large populations of neurons to be broken into cell types. Cell types are either known from the experimental literature or assigned based on cell features like morphology and how they connect to other cell types. Cell typing thus far relied on a mix of varying approaches and generally requires extensive manual intervention by expert annotators. Here, we present a new connectivity-alignment based method for comparing connectomes between hemispheres, individuals and sexes. In contrast to conventional approaches, which heavily rely on morphological comparisons, our method relies solely on topology and is fully automated.  We deploy this method to analyze \~13,000 neurons intrinsic to the ventral nerve cord (VNC) of three adult *Drosophila melanogaster* connectomes (2 male, 1 female). We additionally propose a connectivity-based, unbiased method for defining sex-specific and dimorphic cell types that handles cross-individual variability and reconstruction errors. We discover 203 male-specific, 197 female-specific and 143 sexually dimorphic cell types intrinsic to the ventral nerve cord, the analogue of the vertebrate spinal cord, which integrates peripheral sensory information and drives motor behaviors and other peripheral neural processes. Our analysis yields the first comprehensive census of homologous, sexually dimorphic, and sex-specific neurons intrinsic to the *Drosophila* ventral nerve cord, offering new methodology for more automated comparative connectomics and insights into sex differences in circuits that connect the brain and body.

# **Introduction**

Recent advances in large-scale electron microscopy, machine learning-based automated image segmentation, and technologies for proofreading neural reconstructions have collectively made possible the generation of multiple complete synapse-resolution wiring diagrams, or connectomes, for the *Drosophila melanogaster* model system. These include a whole-brain connectome (Dorkenwald et al. 2024 and Schlegel et al. 2024), a ventral nerve cord connectome [(Takemura et al. 2023\)](https://paperpile.com/c/7p9NJj/E9X8), and, more recently, central nervous system connectomes [(Bates et al. 2025; Berg et al. 2025\)](https://paperpile.com/c/7p9NJj/ETyZ+aKLi). These comprehensive resources open the door to comparative connectomics, systematic comparisons of connectivity across individuals, to uncover biologically meaningful differences. Understanding how circuit architecture varies between individuals, sexes, and species, can reveal which neuronal features are conserved and which are specialized, offering insight into the principles of nervous system organization, the neural basis of behavioral diversity, and the role of evolution in shaping these circuits. For example, discovering that the same circuit is preserved across individuals would suggest strong developmental and functional constraints, indicating that the computation the circuit performs may be more general. In contrast, finding that a circuit differs (through changes in connectivity strength, neuron number, or synaptic partner identity) could reveal mechanisms of variation and plasticity (when comparing individuals), sexual dimorphism (when comparing sexes), or adaptive specialization (when comparing different species).

A key component in comparative connectomics (Costa Nature Methods commentary 2025\) is the ability to accurately align neurons and circuits across different datasets. Neurons within a connectome can be grouped into cell types based on developmental birth order, morphological similarity, and connectivity [(Bates et al. 2019\)](https://paperpile.com/c/7p9NJj/MX9w)(ref). Such grouping is critical not only for compressing a dataset into its functional units to facilitate analysis (for example, the *Drosophila* whole-brain connectome contains 140,000 neurons which can be grouped into just \~10,000 cell types (Dorkenwald et al. 2024\)[(Schlegel et al. 2024\)](https://paperpile.com/c/7p9NJj/sDdt)), but also for making comparisons between datasets. The challenge is establishing homology: reliably determining which individual neuron in a new dataset corresponds to an established cell type in a reference dataset. Well-defined cCell types should largely be conserved across individuals, but the fine-scale morphology and connectivity of neurons within a cell type can vary subtly, making it exceedingly difficult to find the true, topologically equivalent partner across different connectomes without extensive manual expert curation. This problem can be compounded by reconstruction errors in the data (due to manual proofreading (Dorkenwald et al. Nature Methods 2022)). For example, comparing the *Drosophila* hemibrain connectome (\~20K neurons from one female fly (Scheffer et al. 2020)) with the *Drosophila* whole-brain connectome from a second female fly (Zheng et al. 2018\) relied on semi-manual cell type alignment (Schlegel et al. 2024). Once aligned, it was possible to assess the extent of inter-individuality variability in connectivity between the same cell types in the two different female brains. For morphological similarity, NBLAST [(Costa et al. 2016\)](https://paperpile.com/c/7p9NJj/15Ad) serves as the gold standard, but it alone cannot accurately align cell types without human review (refs). Following initial morphology comparisons, connectivity is then used to refine cell types (Scheffer et al. 2020), but the choice of a similarity measure for connectivity has not yet been established. Combining both morphology matching and connectivity comparisons can be laborious, requiring rounds of manual correction, and thus making cell type annotation one of the major bottlenecks in connectomics. 

The main advance of our approach is the complete shift away from using neuronal morphology as the initial scaffolding for alignment. Morphology-first methods inherently face hurdles from noise introduced by non-rigid spatial registration, subtle variation in arbor shape due to both reconstruction error and biological variation, and the need for subjective manual correction to establish confident cell-type homologies. By isolating the alignment process from physical space and relying purely on the network's topological structure, our method establishes a more objective, reproducible, and scalable foundation for comparative analysis. A crucial underlying assumption of a connectivity-first approach is that the aligned parts must be reconstructed with sufficient connectivity data, which is not an absolute requirement for morphology-based mapping methods that do not rely on the surrounding network context. 

Our approach relies on the [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) network alignment method \[ACDC ref\], and we apply it here to the problem of aligning three different connectomic datasets, ventral nerve cords of either male (MANC (ref) or maleCNS (ref)) or female (BANC (ref)). The [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) method alternates between two complementary ways of searching for a good alignment between weighted, directed graphs. One step performs discrete combinatorial search (e.g. swapping pairs of nodes) to increase the alignment score. The other step performs continuous optimization, where the discrete matching problem is relaxed to allow gradient-based methods (such as Frank–Wolfe updates (ref)). By repeatedly alternating between these two steps, the continuous updates help escape local optima reached by the discrete search, while the discrete updates refine the solution to produce a valid permutation that maximizes the final matching score.

Within each dataset we also use the same method to align left and right hemispheres. Comparisons between the left and right hemisphere of each dataset provide a baseline for intra-individual variability. Comparisons between two different male VNC datasets provide a baseline for inter-individual variability.  [![][image2]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) was developed as a result of a FlyWire connectome data challenge (website link). Unlike current methods that rely on subjective, morphology-based constraints, non-rigid spatial registration, and extensive manual expert curation to resolve ambiguous matches and define cell-type boundaries, the proposed method offers both a reduction of bias (by using mathematical calculation of topological variance with simple discrete thresholds to classify dimorphic and sex-specific cell types) and a reduction in human labor (by removing the need for a morphology-first scaffold, pre-identified homologous panels, and labor-intensive human review). [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0), therefore, presents a major step toward developing a more automated, reproducible pipeline for comparative connectomics in future projects. Note however that the [![][image2]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) method, in its current implementation, has a practical limit of approximately 20,000–25,000 neurons. This constraint arises because larger connectivity matrices exceed the memory capacity of standard PCs, making execution time unfeasible. While the VNC intrinsic neurons along with AN/DN neck neurons fit within this limit, a full connectome alignment remains challenging until a more efficient implementation becomes available.

Our application of the [![][image2]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) method automates the discovery of sexually dimorphic and sex-specific neurons (when comparing male and female datasets). Sex differences in the *Drosophila* nervous system are well-documented and known to underlie sex-specific behaviors such as courtship, aggression and reproductive behaviors (refs). These differences arise through the effects of the transcription factors, fruitless and doublesex \- explain. Recent work has identified Fruitless and Doublesex neurons in female and male brain connectomes (Deutsch et al. and Berg et al.) through both EM to LM and EM to EM comparisons. However, we know much less about the sex different circuitry of the ventral nerve cord. Several male-specific cell types of the wing tectulum (important for the production of male courtship song) have been identified through careful EM/LM matching, and their circuitry analyzed within the MANC connectome (Lilvis et al. Current Biology), but to date there has been no comprehensive comparison between the connectomes of male and female ventral nerve cords. This presents a major gap in our understanding of the differences in circuitry between male and female nervous systems that connect the brain and body. 

The adult whole-brain connectome recovered 68% of types initially defined in the partial hemibrain connectome[6](https://www.nature.com/articles/s41592-025-02946-2#ref-CR6),[7](https://www.nature.com/articles/s41592-025-02946-2#ref-CR7),[8](https://www.nature.com/articles/s41592-025-02946-2#ref-CR8). The recent whole-CNS fly connectome has so far cross-matched 74% of neurons to existing connectomes[11](https://www.nature.com/articles/s41592-025-02946-2#ref-CR11). Comparisons across the sexes in *C. elegans* and in the fly (though with limited data so far) have shown high conservation of cell types, with no more than 5% of neurons being sex specific and 5% dimorphic.

# **Results**

**Comparison with prior work**: Current comparative connectomics pipelines typically operate through a phased, morphology-first framework, with the assumption that morphological, spatial, or transcriptomic information provides the primary scaffold for cell identity (refs). There are only two species (so far) for which there exists multiple connectome datasets, and for which comparisons are possible: *C. elegans* and *D. melanogaster*. In *C. elegans*, cell types are largely identified by the invariant position of the soma (Cook et al. 2025). In *Drosophila*, cell body position is not sufficiently stereotyped, and the morphology of neuronal arbors is used (ref). Matching cells by morphology typically involves non-rigid spatial transformations (registering the datasets to a common coordinate space), to facilitate comparisons with existing connectomes (Schlegel et al. 2024\) or LM (light microscopy) data (Deutsch et al. 2025 and Sturner et al. 2025). Then, highly confident morphological matches are identified utilizing pairwise NBLAST scores coupled with the Hungarian (linear assignment) algorithm and co-visualization with human verification. Because the best match for each cell is often unclear (due to many equally good matches with similar NBLAST scores), this process is often executed in iterative rounds with progressively decreasing confidence thresholds, and its established seed of confident matches subsequently serves as a shared coordinate system. Connectivity vectors are then constructed based on the identity of synaptic partners within this confident set, and cosine similarity is employed to quantify the resulting "connectivity distance" between neurons. Finally, the defined connectivity distance is used to resolve the remaining unmatched neurons with ambiguous / low confidence morphology-based matches. This annotation process can be slow and represents a major bottleneck in the generation of connectomes.

Our proposed method diverges from prior work in three ways:

1\. **Connectivity based initial alignment**: We dispense with the use of morphology, and matching is executed purely on synaptic connectivity utilizing the [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) graph alignment algorithm. This isolates the matching process from the physical space, entirely bypassing the potential noise, registration errors, and measurement gaps inherent to non-rigid spatial transformations and NBLAST score calculations.   
2\. **Continuous alignment weighting**: When building connectivity features, we do not rely on a discretely bounded seed of "confident" partner identities. Instead, we systematically incorporate the alignment score of every matched pair into the feature space. Feature coordinates are defined by multiplying the connection strength (normalized synapse count) to a partner by that partner's corresponding [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) alignment score, allowing the entire network topology to dynamically weight the mapping.  
3\. **Weighted Intersection as the similarity metric**: To evaluate the similarity between the resulting connectivity feature vectors, the proposed method utilizes weighted Jaccard similarity  / weighted intersection (Matsliah et al., 2024). This metric is explicitly selected over cosine similarity (which discards critical magnitude information) and inverse L2 distance (which is disproportionately dominated by the largest values). This choice provides a more robust and accurate quantification of topological equivalence in high-dimensional, sparse connectomic data, as evident in experiments summarized in Section XXX.  

Our approach is grounded in the theoretical premise that synaptic connectivity alone is sufficient for defining neuronal cell types, and that morphological features are eventually a proxy for connectivity (see discussion in NTAC paper XXX).

![][image3]  
**Fig. 1 | Connectome alignment and weighted feature embedding** Visual depiction of the proposed connectivity-alignment-first pipeline, illustrating the process of mapping neuronal adjacency matrices and constructing dynamically weighted feature vectors to evaluate similarity. **(a)** Network level view of the preprocessing and alignment workflow. The network pairs depict connectomes A and B before synapse weight normalization (top), after normalization (middle), and after alignment of B to A (bottom). **(b)** Initial alignment of adjacency matrices via the [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) graph alignment algorithm. The sparse adjacency matrix of connectome A is matched to the adjacency matrix of connectome B. [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) optimizes for weighted intersection, outputting a permutation matrix [![][image4]](https://www.codecogs.com/eqnedit.php?latex=%5Csigma#0) that reorders the neurons in B to maximize topological correspondence with A. **(c)** Construction of dynamically weighted connectivity feature vectors for a pair of neurons [![][image5]](https://www.codecogs.com/eqnedit.php?latex=a_x%20%5Cin%20A#0) and [![][image6]](https://www.codecogs.com/eqnedit.php?latex=b_y%20%5Cin%20B#0). Each coordinate represents the number of input or output synapses (horizontal bars) to one of upstream or downstream partners respectively. The connectivity similarity between a query neuron [![][image7]](https://www.codecogs.com/eqnedit.php?latex=a_x#0) a candidate neuron [![][image8]](https://www.codecogs.com/eqnedit.php?latex=b_y#0) scales each coordinate by the [![][image1]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) alignment score [![][image9]](https://www.codecogs.com/eqnedit.php?latex=%5Comega_i#0). This mechanism adjusts the influence of each synaptic connection during similarity evaluation. Note: while in this example the number of nodes match exactly, in real applications any imbalances are corrected by adding singleton (disconnected) padding nodes.

**Identifying sex-different circuitry between male and female datasets**: In (Berg et al.) and (Stürner et al.), cell types were categorized into isomorphic, sexually dimorphic, and sex-specific, based on structural heuristics and manual expert curation. Sex-specific neurons were defined as groups exhibiting bilateral consistency within the same sex but lacking a probable morphological match in the opposite-sex dataset. Sexually dimorphic neurons were defined as having a shared developmental origin and cross-sex matches, but exhibiting consistent morphological differences, such as extra branches, that resulted in corresponding changes in connectivity. In (Cachero et al.), sexual dimorphism was evaluated by mapping sex-specific transcriptional profiles to the connectome. This approach identified female-specific apoptosis and transcriptional divergence as the primary global drivers of sex-specification.

**In contrast, our method defines sex differences quantitatively.** For each neuron, we compute the ratio between (i) similarity to its cross-hemisphere twin and (ii) similarity to its closest cross-dataset match. This cross-hemisphere/cross-dataset similarity ratio provides a continuous measure of dimorphism (Fig. XXX).

To assign discrete labels, we threshold the median ratio within each cell type:

* **\< 2** → isomorphic

* **2–4** → sexually dimorphic

* **\> 4** → sex-specific

These thresholds were calibrated using previously characterized cell types.

By deriving dimorphism status directly from the weighted feature space of the graph alignment, this approach replaces visual morphological evaluation, transcriptomic profiling, and human judgment with discrete mathematical ratios. This provides a streamlined analytical pipeline that minimizes the potential for subjective bias in classification.

While our method is piloted on the intrinsic neurons of the *Drosophila* ventral nerve cord (\~13,000 neurons), it can in principle be applied to any proofread connectome dataset for which a reference dataset is available. 

![][image10]  
**Fig. 2 | Classification and spatial distribution of sex specific and dimorphic neurons** We utilized the connectivity-based similarity ratio (cross-hemisphere twin similarity divided by cross-dataset closest match similarity) to systematically categorize VNC cell types into isomorphic, sexually dimorphic, and sex-specific populations. **(a)** Dimorphism ratios for known cell types (MANC to BANC): Scatter plot evaluating the similarity ratios of male neurons mapped to the female connectome. This panel highlights cell types with previously established classifications to validate the thresholds against biological ground truth. **(b)** Dimorphism ratios for known cell types (BANC to MANC): The inverse mapping of female neurons to the male connectome, evaluating previously classified cell types. **(c, d)** Classification of dimorphic and sex specific cell types with the established thresholds. **(e, f)** Spatial heatmaps of male (e) and female (f) dimorphisms: Visualization of neuron centroids colored by their computed dimorphism ratio to highlight highly dimorphic and sex-specific regions. **(g)** Intra-sex baseline variability (MANC to Male CNS): Scatter plot of similarity ratios comparing two independent male nerve cord datasets. This male-to-male comparison establishes a quantitative baseline for inter-individual network variation. **(h)** Cross-dataset dimorphism validation (BANC to Male CNS): Spatial heatmap of female neuron centroids colored by their dimorphism ratio when mapped against a second, independent male dataset, confirming the structural localizations observed in the primary comparison.

![][image11]

**Fig. 3 | Stereotypy and variation of cell populations and synaptic connectivity across hemispheres and datasets** Evaluation of the network stereotypy of VNC-intrinsic cell types by comparing total cell counts and input-normalized synaptic fractions cross hemispheres and datasets. **(a–c)** Comparing the number of cells per type in the left versus the right hemisphere for female (BANC) and male (MANC, Male CNS) datasets. High correlation indicates strong bilateral symmetry in population sizes. **(d-f)** Comparison of total cell counts per type across datasets. Significant deviations from the diagonal highlight numerically dimorphic or sex-specific cell populations. **(g–i)** Scatter plots (log scale) comparing the input-normalized fraction of synapses between cell types in the left versus the right hemisphere for BANC, MANC and Male CNS. Strong correlation along the diagonal demonstrates bilateral symmetry in synaptic weights. **(j-l)** Comparison of type-to-type synaptic fractions across sexes. Connections are stratified based on the dimorphism categories of the participating cell types.

![][image12]

**Fig. 5 | Synaptic connectivity bias of sexually dimorphic types** Each point represents a distinct VNC intrinsic cell type, plotted by the percentage of its input/output synapses (x-axis) and connections (y-axis) with sex specific cell types. Points are colored by the cell type's own classification: isomorphic (grey), dimorphic (green), or sex-specific (red). Labeled points highlight cell types with exceptional biases. **(a)** BANC VNC intrinsic cell types upstream connectivity with female specific cells. Isomorphic types cluster near the origin, while dimorphic and sex-specific types exhibit significantly higher connectivity with female specific partners. **(b)** same but now with respect to downstream partners. **(c, d)** A parallel analysis in the male MANC VNC connectome.

![][image13]  
**Fig. 6** | **Distribution of sex specific and dimorphic neurons in developmental hemilineages** The stacked bar charts display dimorphic and sex-specific cells (top segment) together with isomorphic cells (base segment) for each hemilineage. Neurons with \>50% of their synapses in the abdominal ganglion are grouped in ABg and the vast majority of them lack hemilineage assignment in both datasets. Remaining neurons without hemilineage assignment are grouped in XXX.Text annotations denote the absolute number of dimorphic cells, followed by the percentage of the total hemilineage population that exhibits dimorphism. **(a)** Dimorphism distribution in the male VNC (MANC). **(b)** Dimorphism distribution in the female VNC (BANC).

# **Methods**

## **Datasets**

The comparison between the MaleCNS VNC and MANC datasets is foundational for establishing a reliable baseline of inter-individual variability in the Drosophila ventral nerve cord. Both connectomes originate from different individuals of the same male genotype and were reconstructed using the same technique (FIB-SEM imaging, followed by automated synapse detection and manual proofreading, \[XXX\]). By quantifying the topological differences between these two individuals, we establish an expected level of natural variation that is independent of sex. This baseline is critical, as it allows us to confidently classify the observed topological differences between a male (e.g., MANC) and a female (BANC, which was acquired with TEM) as a true sexual dimorphism only if the difference significantly exceeds the established inter-individual variation threshold.

**Fig. 7 | Anatomy and annotation coverage of the Drosophila central nervous system connectomes** Overview of the gross anatomy of the Drosophila central nervous system and the quantitative distribution of reconstructed neurons in the VNC. **(a)** Anatomical schematic of the Drosophila central nervous system (CNS). The diagram illustrates the major structural divisions, highlighting the spatial relationship between the brain and the ventral nerve cord (VNC), which are physically linked by the neck connective. **(b)** Distribution of reconstructed neurons and cell types across VNC connectomes. The grouped bar chart quantifies the total number of reconstructed neurons across four primary functional categories (Neck, VNC Sensory, VNC Motor, and VNC Intrinsic) for the female connectome (BANC), the primary male connectome (MANC), and the secondary male connectome (Male CNS). The absolute cell count is indicated at the apex of each bar, with the corresponding number of unique cell types assigned to each population annotated vertically. This establishes the biological scale, composition, and typing completeness of the datasets evaluated in the comparative alignment pipeline.

## **Alignment and Similarity Score definitions**

Alignment weight [![][image14]](https://www.codecogs.com/eqnedit.php?latex=%5Comega_i#0) is defined as the weighted Jaccard similarity between the [![][image15]](https://www.codecogs.com/eqnedit.php?latex=%5Csigma#0)\-permuted in- and out-synaptic neighborhoods of matched node pairs [![][image16]](https://www.codecogs.com/eqnedit.php?latex=a_i#0) and [![][image17]](https://www.codecogs.com/eqnedit.php?latex=b_%7B%5Csigma\(i\)%7D#0).  
The similarity score between any pair of neurons [![][image18]](https://www.codecogs.com/eqnedit.php?latex=a_x#0) and [![][image19]](https://www.codecogs.com/eqnedit.php?latex=b_y#0) is then computed as a  weighted Jaccard similarity of their directional synaptic feature vectors, scaled by coordinate-wise alignment [![][image20]](https://www.codecogs.com/eqnedit.php?latex=%5Comega_i#0) of each partner.

Alignment score for [![][image21]](https://www.codecogs.com/eqnedit.php?latex=%5Csigma#0) \- defined for [![][image22]](https://www.codecogs.com/eqnedit.php?latex=N#0) matched pairs of nodes [![][image23]](https://www.codecogs.com/eqnedit.php?latex=a_i#0), [![][image24]](https://www.codecogs.com/eqnedit.php?latex=b_%7B%5Csigma\(i\)%7D#0)  
[![][image25]](https://www.codecogs.com/eqnedit.php?latex=%5Comega_i%20%3D%20%5Cfrac%7B%5Csum_%7Bj%3D1%7D%5EN%20%5CBig%5B%5Cmin\(%5Cmathrm%7Bsyn%7D\(a_j%2Ca_i\)%2C%5Cmathrm%7Bsyn%7D\(b_%7B%5Csigma\(j\)%7D%2Cb_%7B%5Csigma\(i\)%7D\)\)%2B%5Cmin\(%5Cmathrm%7Bsyn%7D\(a_i%2Ca_j\)%2C%5Cmathrm%7Bsyn%7D\(b_%7B%5Csigma\(i\)%7D%2Cb_%7B%5Csigma\(j\)%7D\)\)%5CBig%5D%7D%7B%5Csum_%7Bj%3D1%7D%5EN%20%5CBig%5B%5Cmax\(%5Cmathrm%7Bsyn%7D\(a_j%2Ca_i\)%2C%5Cmathrm%7Bsyn%7D\(b_%7B%5Csigma\(j\)%7D%2Cb_%7B%5Csigma\(i\)%7D\)\)%2B%5Cmax\(%5Cmathrm%7Bsyn%7D\(a_i%2Ca_j\)%2C%5Cmathrm%7Bsyn%7D\(b_%7B%5Csigma\(i\)%7D%2Cb_%7B%5Csigma\(j\)%7D\)\)%5CBig%5D%7D#0)

Similarity score for [![][image26]](https://www.codecogs.com/eqnedit.php?latex=%5Csigma#0) \- defined for any pair of nodes [![][image27]](https://www.codecogs.com/eqnedit.php?latex=a_x#0), [![][image28]](https://www.codecogs.com/eqnedit.php?latex=b_y#0)  
[![][image29]](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BSim%7D\(a_x%2Cb_y\)%3D%5Cfrac%7B%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%5Comega_i%20%5CBig%5B%5Cmin\(a_i%5E%7B%5Cmathrm%7Bin%7D%7D%2Cb_%7B%5Csigma\(i\)%7D%5E%7B%5Cmathrm%7Bin%7D%7D\)%2B%5Cmin\(a_i%5E%7B%5Cmathrm%7Bout%7D%7D%2Cb_%7B%5Csigma\(i\)%7D%5E%7B%5Cmathrm%7Bout%7D%7D\)%5CBig%5D%7D%7B%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20%5Comega_i%20%5CBig%5B%5Cmax\(a_i%5E%7B%5Cmathrm%7Bin%7D%7D%2Cb_%7B%5Csigma\(i\)%7D%5E%7B%5Cmathrm%7Bin%7D%7D\)%2B%5Cmax\(a_i%5E%7B%5Cmathrm%7Bout%7D%7D%2Cb_%7B%5Csigma\(i\)%7D%5E%7B%5Cmathrm%7Bout%7D%7D\)%5CBig%5D%7D#0)

where  
[![][image30]](https://www.codecogs.com/eqnedit.php?latex=a_i%5E%7B%5Cmathrm%7Bin%7D%7D%20%3A%3D%20%5Cmathrm%7Bsyn%7D\(a_i%2Ca_x\)%2C%20%5Cquad%20a_i%5E%7B%5Cmathrm%7Bout%7D%7D%20%3A%3D%20%5Cmathrm%7Bsyn%7D\(a_x%2Ca_i\)#0)  
[![][image31]](https://www.codecogs.com/eqnedit.php?latex=b_i%5E%7B%5Cmathrm%7Bin%7D%7D%20%3A%3D%20%5Cmathrm%7Bsyn%7D\(b_i%2Cb_y\)%2C%20%20%5Cquad%20b_i%5E%7B%5Cmathrm%7Bout%7D%7D%20%3A%3D%20%5Cmathrm%7Bsyn%7D\(b_y%2Cb_i\)#0)

To identify dimorphisms, we compare two quantities for each neuron: its similarity to its closest cross-hemisphere twin (within the same dataset) versus its similarity to its closest cross-dataset match. The median ratio of these similarities is computed as the final similarity ratio for the cell type. This ratio provides a quantitative measure of topological variance, where a value of 1 indicates high homology. Neuron types are then classified using two discrete thresholds: a type is classified as sexually dimorphic if its median similarity ratio is between 2 and 4, and as sex-specific if the ratio exceeds 4\. These thresholds were not arbitrarily chosen but were empirically set to minimize the false-positive classification rate between the two male datasets (MaleCNS VNC and MANC), thereby establishing a conservative baseline for defining true sex differences.

##  **Effectiveness of Similarity Measures** 

To evaluate effectiveness of similarity measures we considered: [![][image32]](https://www.codecogs.com/eqnedit.php?latex=Cosine#0) similarity, Inverse [![][image33]](https://www.codecogs.com/eqnedit.php?latex=L_1#0) distance, Inverse [![][image34]](https://www.codecogs.com/eqnedit.php?latex=L_2#0) distance, Weighted Intersection and Weighted Jaccard similarity. In all experiments, we generated sparse vectors for each neuron to represent its connectivity profile, including both its input and output partners (or partner types) and the corresponding synaptic weights. For connectivity vectors [![][image35]](https://www.codecogs.com/eqnedit.php?latex=u#0) and [![][image36]](https://www.codecogs.com/eqnedit.php?latex=v#0), the evaluated similarity measures are defined as follows:

[**![][image37]**](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BCosine%7D\(u%2Cv\)%20%3D%20%5Cfrac%7Bu%20%5Ccdot%20v%7D%7B%5C%7Cu%5C%7C%5C%2C%5C%7Cv%5C%7C%7D#0)  
[![][image38]](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BInverse%5C%20L1%7D\(u%2Cv\)%20%3D%20%5Cfrac%7B1%7D%7B1%20%2B%20%5Csum_i%20%7Cu_i%20-%20v_i%7C%7D#0)  
[![][image39]](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BInverse%5C%20L2%7D\(u%2Cv\)%20%3D%5Cfrac%7B1%7D%7B1%20%2B%20%5Csqrt%7B%5Csum_i%20\(u_i%20-%20v_i\)%5E2%7D%7D#0)  
[![][image40]](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BWeighted%5C%20Intersection%7D\(u%2Cv\)%20%3D%5Csum_i%20%5Cmin\(u_i%2C%20v_i\)#0)  
[![][image41]](https://www.codecogs.com/eqnedit.php?latex=%5Cmathrm%7BWeighted%5C%20Jaccard%7D\(u%2Cv\)%20%3D%5Cfrac%7B%5Csum_i%20%5Cmin\(u_i%2C%20v_i\)%7D%7B%5Csum_i%20%5Cmax\(u_i%2C%20v_i\)%7D#0)

**Experiment 1: Intra-Type Retrieval Task:** Our first evaluation was a retrieval task to test how well a measure can identify neurons of the same type. For each dataset, we conducted 1,000 trials, selecting a random neuron and finding its most similar partner using each of the four measures. We recorded success when the top hit matched the same neuronal type.

**Experiment 2: Cross-Hemisphere Sub-circuit Alignment Task:** Our second evaluation directly tested the ability of each measure to find the identity (*exactly* *correct)* alignment (graph matching) between homologous sub-circuits. We sampled 1,000 connected 10-neuron subcircuits from the left hemisphere and their ground-truth counterparts on the right. For each measure, we evaluated all *10\!* possible permutations and recorded a success if the measure's score was maximized at the identity permutation.

**Experiment 3: Large-Scale Alignment via Data Challenge:** After selecting Weighted Jaccard as the best performing similarity measure, we hosted a public data challenge (XXX ref) attracting dozens of participants including the creators of state of the art network alignment methods (FAQ \[XXX\], SANA \[XXX\]). The challenge ran from October 2024 to February 2025 \[XXX\], and the winning algorithm [![][image42]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) achieved top-scoring alignment between BANC VNC and MANC.

The winning [![][image43]](https://www.codecogs.com/eqnedit.php?latex=AC%20%5Coplus%20DC#0) algorithm, using Weighted Intersection (which is the simpler equivalent to Weighted Jaccard for this task \- see proof in Appendix XXX) as its objective function, was then verified on large scale alignment tasks, comparing FlyWire, BANC, Male CNS, MANC and Hemibrain datasets. Quality of alignment is measured by counting the fraction of aligned pairs that share the same type across datasets. 

**Fig. 8 | Evaluating similarity measures** We assessed the accuracy and effectiveness of different similarity metrics for neuron connectivity. **(a)** Intra-dataset neuron retrieval: For each dataset and measure, we performed 1,000 trials. In each trial, a random query neuron (excluding untyped and singleton cells) was matched to its closest neighbor using feature vectors composed of input and output synapse counts grouped by partner types. Error rates are reported as the percentage of trials where the closest match was not the same cell type as the query. **(b)** Evaluating similarity measures for cross-hemisphere sub-circuit alignment: To test the ability of different metrics to identify the exact topological correspondence between homologous networks, we sampled 1,000 pairs of isomorphic, connected 10-neuron sub-circuits across the left and right hemispheres. To ensure an unambiguous ground truth, these sub-circuits were composed exclusively of singleton cell types (neurons with exactly one instance per hemisphere). For each sub-circuit pair, we evaluated all possible mappings ($10\!$ permutations). We compared Weighted Jaccard and Cosine similarity; as demonstrated in the Methods, for the specific mathematical task of optimizing permutations, Weighted Jaccard is equivalent to inverse L1 distance and Weighted Intersection, while Cosine similarity is equivalent to inverse L2 distance. Error rates represent the percentage of the 1,000 trials where the similarity metric failed to be globally maximized by the known ground-truth alignment. Lower error percentages indicate superior alignment accuracy.

## **Equivalence of Similarity Measures for Network Alignment**

For the specific task of finding the best alignment (i.e., finding the optimal permutation [![][image44]](https://www.codecogs.com/eqnedit.php?latex=P#0) of a neuron's partners), two of the considered measure pairs become equivalent.

Equivalence of Cosine Similarity and Inverse L2 Distance

The cosine similarity for permutation [![][image45]](https://www.codecogs.com/eqnedit.php?latex=P#0) is

[![][image46]](https://www.codecogs.com/eqnedit.php?latex=%5Ccos\(P\)%20%3D%20%5Cfrac%7B%5Clangle%20x_P%2C%20y%20%5Crangle%7D%7B%5C%7Cx_P%5C%7C_2%20%5C%2C%20%5C%7Cy%5C%7C_2%7D.#0)

Because [![][image47]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%5C%7C_2%20%3D%20%5C%7Cx%5C%7C_2#0) for all [![][image48]](https://www.codecogs.com/eqnedit.php?latex=P#0), maximizing [![][image49]](https://www.codecogs.com/eqnedit.php?latex=%5Ccos\(P\)#0) is equivalent to maximizing the inner product

[![][image50]](https://www.codecogs.com/eqnedit.php?latex=%5Clangle%20x_P%2C%20y%20%5Crangle.#0)

The squared Euclidean distance is

[![][image51]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%20-%20y%5C%7C_2%5E2%20%3D%20%5C%7Cx_P%5C%7C_2%5E2%20%2B%20%5C%7Cy%5C%7C_2%5E2%20-%202%5Clangle%20x_P%2C%20y%20%5Crangle.#0)

The first two terms are constant with respect to [![][image52]](https://www.codecogs.com/eqnedit.php?latex=P#0), so minimizing [![][image53]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%20-%20y%5C%7C_2#0) (or maximizing its inverse) is also equivalent to maximizing [![][image54]](https://www.codecogs.com/eqnedit.php?latex=%5Clangle%20x_P%2C%20y%20%5Crangle#0).

Therefore, the permutation that maximizes cosine similarity is exactly the same as the one that minimizes [![][image55]](https://www.codecogs.com/eqnedit.php?latex=L_2#0) distance (or maximizes inverse [![][image56]](https://www.codecogs.com/eqnedit.php?latex=L_2#0)).

## **Equivalence of Inverse L1 Distance and Weighted Intersection**

The [![][image57]](https://www.codecogs.com/eqnedit.php?latex=L_1#0) distance for permutation [![][image58]](https://www.codecogs.com/eqnedit.php?latex=P#0) is

[![][image59]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%20-%20y%5C%7C_1%20%3D%20%5Csum_i%20%7Cx_%7BP%2Ci%7D%20-%20y_i%7C.#0)

For any real numbers [![][image60]](https://www.codecogs.com/eqnedit.php?latex=u%2Cv#0),

[![][image61]](https://www.codecogs.com/eqnedit.php?latex=%7Cu%20-%20v%7C%20%3D%20u%20%2B%20v%20-%202%5Cmin\(u%2C%20v\).#0)

Applying this identity componentwise gives

[![][image62]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%20-%20y%5C%7C_1%20%3D%20%5Csum_i%20x_%7BP%2Ci%7D%20%2B%20%5Csum_i%20y_i%20-%202%20%5Csum_i%20%5Cmin\(x_%7BP%2Ci%7D%2C%20y_i\).#0)

The first two sums are independent of [![][image63]](https://www.codecogs.com/eqnedit.php?latex=P#0), so minimizing [![][image64]](https://www.codecogs.com/eqnedit.php?latex=%5C%7Cx_P%20-%20y%5C%7C_1#0) (or equivalently maximizing its inverse) is equivalent to maximizing

[![][image65]](https://www.codecogs.com/eqnedit.php?latex=%5Csum_i%20%5Cmin\(x_%7BP%2Ci%7D%2C%20y_i\)%2C#0)

the *weighted intersection* (or coordinate-wise overlap).

Hence, the permutation that minimizes [![][image66]](https://www.codecogs.com/eqnedit.php?latex=L_1#0) distance is identical to the one that maximizes the weighted intersection score.

## **Equivalence of Weighted Intersection and Weighted Jaccard**

Let [![][image67]](https://www.codecogs.com/eqnedit.php?latex=A%3D\(a_%7Bij%7D\)#0) and [![][image68]](https://www.codecogs.com/eqnedit.php?latex=B%3D\(b_%7Bij%7D\)#0) be weighted directed adjacency matrices with  
[![][image69]](https://www.codecogs.com/eqnedit.php?latex=a_%7Bij%7D%2C%20b_%7Bij%7D%20%5Cge%200#0). Let [![][image70]](https://www.codecogs.com/eqnedit.php?latex=P#0) be a permutation representing a node alignment,  
and define the permuted matrix [![][image71]](https://www.codecogs.com/eqnedit.php?latex=b%5E%7BP%7D_%7Bij%7D%20%3D%20b_%7BP\(i\)P\(j\)%7D#0).

Define the weighted intersection  
[![][image72]](https://www.codecogs.com/eqnedit.php?latex=I\(P\)%20%3D%20%5Csum_%7Bi%2Cj%7D%20%5Cmin\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)#0)  
and the weighted Jaccard similarity  
[![][image73]](https://www.codecogs.com/eqnedit.php?latex=J\(P\)%20%3D%5Cfrac%7B%5Csum_%7Bi%2Cj%7D%20%5Cmin\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)%7D%7B%5Csum_%7Bi%2Cj%7D%20%5Cmax\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)%7D#0)

For any nonnegative [![][image74]](https://www.codecogs.com/eqnedit.php?latex=x%2Cy#0) we have  
[![][image75]](https://www.codecogs.com/eqnedit.php?latex=%5Cmax\(x%2Cy\)%20%3D%20x%20%2B%20y%20-%20%5Cmin\(x%2Cy\)#0)

Applying this identity elementwise and summing gives

[![][image76]](https://www.codecogs.com/eqnedit.php?latex=%5Csum_%7Bi%2Cj%7D%20%5Cmax\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)%3D%5Csum_%7Bi%2Cj%7D%20a_%7Bij%7D%2B%5Csum_%7Bi%2Cj%7D%20b%5E%7BP%7D_%7Bij%7D-%5Csum_%7Bi%2Cj%7D%20%5Cmin\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)#0)

Let

[![][image77]](https://www.codecogs.com/eqnedit.php?latex=S_A%20%3D%20%5Csum_%7Bi%2Cj%7D%20a_%7Bij%7D%2C%20%5Cqquad%20S_B%20%3D%20%5Csum_%7Bi%2Cj%7D%20b_%7Bij%7D#0)

Since permutation preserves sums,

[![][image78]](https://www.codecogs.com/eqnedit.php?latex=%5Csum_%7Bi%2Cj%7D%20b%5E%7BP%7D_%7Bij%7D%20%3D%20S_B#0)

Therefore  
[![][image79]](https://www.codecogs.com/eqnedit.php?latex=%5Csum_%7Bi%2Cj%7D%20%5Cmax\(a_%7Bij%7D%2C%20b%5E%7BP%7D_%7Bij%7D\)%3DS_A%20%2B%20S_B%20-%20I\(P\)#0).

Substituting into the definition of [![][image80]](https://www.codecogs.com/eqnedit.php?latex=J\(P\)#0) yields  
[![][image81]](https://www.codecogs.com/eqnedit.php?latex=J\(P\)%20%3D%20%5Cfrac%7BI\(P\)%7D%7BS_A%20%2B%20S_B%20-%20I\(P\)%7D#0)

Let [![][image82]](https://www.codecogs.com/eqnedit.php?latex=C%20%3D%20S_A%20%2B%20S_B#0), which is constant with respect to [![][image83]](https://www.codecogs.com/eqnedit.php?latex=P#0). Then

[![][image84]](https://www.codecogs.com/eqnedit.php?latex=J\(P\)%20%3D%20f\(I\(P\)\)%2C%20%5Cqquad%20f\(I\)%20%3D%20%5Cfrac%7BI%7D%7BC-I%7D#0)

The derivative is  
[![][image85]](https://www.codecogs.com/eqnedit.php?latex=f'\(I\)%20%3D%20%5Cfrac%7BC%7D%7B\(C-I\)%5E2%7D%20%3E%200#0)

for [![][image86]](https://www.codecogs.com/eqnedit.php?latex=0%20%5Cle%20I%20%3C%20C#0), so [![][image87]](https://www.codecogs.com/eqnedit.php?latex=f#0) is strictly increasing. Hence for any permutations [![][image88]](https://www.codecogs.com/eqnedit.php?latex=P_1%2C%20P_2#0),

[![][image89]](https://www.codecogs.com/eqnedit.php?latex=I\(P_1\)%20%3E%20I\(P_2\)%5Ciff%20J\(P_1\)%20%3E%20J\(P_2\)#0).

Therefore the two objectives induce the same ordering over permutations and

[![][image90]](https://www.codecogs.com/eqnedit.php?latex=%5Carg%5Cmax_P%20I\(P\)%20%3D%20%5Carg%5Cmax_P%20J\(P\)#0)

Thus maximizing weighted intersection is equivalent to maximizing weighted  
Jaccard similarity for adjacency matrix alignment.

# **Analysis of VNC Sexual Dimorphism**

**![][image91]**  
**Figure 9 (Proposed Figure 3\) |** Sex-Specific and Dimorphic Neurons in the VNC. A) Schematic of cell super classes in the VNC (naming conventions correspond to Bates et al., 2026). B) Comparison of the numbers of cells and cell types in each super class. C) Cells that have been fully or partially reconstructed in MANC and BANC VNC that are as of yet untyped. Cells without a final super class assigned are indicated in by the empty portion of the bars. D) Counts of cells and cell types across all super classes broken down by their dimorphism. E) Counts of INs and IN types by dimorphism. F) Kernel density estimate of synapse density based on postsynapse coordinates presented in JRC2018 reference brain space. Arrowheads indicate the mesothoracic triangle.

Possible addition to Proposed Figure 3 for next itertion, with suggestion from Alex

1. Network layout of sex-specific and dimorphic cell types in BANC and MANC (organized by top inputs and outputs) \- [Arie Matsliah](mailto:am1736@princeton.edu)  
   1. What population?  
   2. Try to tie things together with ANs and DNs  
      1. Need to finish female-specific ANs 

![][image92]  
**Figure 10 (Proposed Figure SX) |** A) Absolute numbers of dimorphisms in all VNC neurons in BANC and MANC, excluding SNs. SNs are excluded because their numbers vary widely between BANC and MANC and matching is still incomplete. B) Absolute numbers of dimorphisms in all VNC INs.

## **Sex Dimorphisms in the Wing Neuropil**

In this section compare male song circuit neurons (MANC) to female \- analyze connectivity within the circuit as well as to DNs, ANs, sensory neurons, and motor neurons. Determine what function the sex dimorphic neurons perform in females. And determine what happens (in females) to isomorphic neurons that receive a lot of song sex-specific input in male.

Figure 4: Circuits of the Wing Tectulum:

Key Observations:

* Analysis of MANC connectome reveals expanded putative song circuit (panels A and B)  
* Sex Dimorphic and Sex Specific neurons tend to not be directly presynaptic to wing motor neurons (panels B, C, D)  
* When comparing the part of the song circuit present in both males and females, these circuit diagrams are still different and this is due to both edges from dimorphic and ismorphic neurons (panel C and D)  
* When comparing flight circuitry (flight steering and flight power), connectivity and circuit diagram is similar in males and females for the same neuropil (wing) \- Fig 4  Supp Fig 1  
* Some isomorphic neurons (like AN02A001) have dimorphic outputs in male and female \- for this particular neuron this is due to an absence of male-specific neurons in the wing neuropil in females. Inhibitory AN can exert larger influence in female wing circuits missing male specific inputs to the same neurons. (Fig 4 Supp Fig 2\)  
* Dimorphic song circuit neurons include (list) \- they have different morphology even though they were identified by connectivity only \- and when comparing their inputs and outputs, we find that they have stronger inhibition in female, and lack female-specific outputs (versus in male, where they have lots of male specific targets). Any new function identified? If so put it here (panel e)  
* There is only one female-specific wing neuropil cell type (INT2996) that is an output of the shared song circuit \- we analyzed its conenctivity in females \- it is inhibitory (Glu) and we find that it functions in the DNa04 flight steering pathway, and its function seems to be to inhibit the hg muscles \- but figure out why\! (Fig 4 \- Supp Fig 3\) \- this figure needs a bit more work  
    
1. Core Song Circuit in MANC \- Neurons identified in Lillvis plus strongly connected neurons  
2. This same circuit in BANC  
3. Compare with flight circuits (from Cheong et al. \- compare BANC and MANC \- fewer sex differences in same neuropil?)  
4. New male-specific and sex-dimorphic song neurons (in MANC, but not studied in Lillvis) \- proposal for what they do functionally? Here analyze connectivity to MNs, SNs, DNs, and ANs to make predictions (and to other song circuit neurons with known functions).  
5. For sex-dimorphic song circuit neurons (present in BANC) \- what do they do in female VNC? Here analyze connectivity to MNs, SNs, DNs, and ANs to make predictions. Are connections largely reduced (so these neurons don’t strongly influence the network) or do they connect with different cell types \- here the focus is on function in the female nerve cord  
6. For male-specific neurons in MANC (like pMP2, pIP10, TN1A, etc.) \- focus on their non sex-dimorphic partners (that get a lot of input from male-specific neurons in MANC) \- what happens to these neurons in BANC? All other connections stay the same (most likely if they get matched as the same cell type) \- is there some compensation for the loss of these connections? Look closely at where the male-specific synapses are and what happens at these locations in the BANC neurons.

Figure 5: Sex Dimorphisms in the Abdominal Ganglion:

Key Observations:

* Identified all female, male specific and dimorphic neurons in the abdominal ganglia (panel a)  
* Stats on how many neurons of each type there are in the AbG versus other neuropils of the VNC (panel b)  
* DNp13 studied in Sturner et al. (compared MANC to FANC) but FANC was missing AbG \- so here we repeat the analysis but with a complete VNC (panel c). Conclusions: ; Contrast with Sturner et al \- they predict a role for DNp13 in male song, but make note that our circuit diagram tells a different story \- DNp13 is only a tiny fraction of TN1A’s inputs \- rather DNp13’s strongest influence in song circuitry is on vms16….in this section emphasize the contrast with Sturner \- complete VNC connectome needed for interpretation and also one should examine the strength of influence of neuron on its output (relative to other inputs that neuron gets)  
* oviDNx figure \- look for shared outputs of oviDNx and other oviDNs

Figure 5: Circuits of the Abdominal Ganglion:  
For this figure, provide stats and overview of the cell types we have labeled as sex-specific or sex-dimorphic (similar to Fig 3, but focus on ANm) \- and provide 1 deep dive (cell types downstream of oviDNs (drives egg laying), DNp13 (drives ovipositor extrusion, a rejection behavior), and vpoDN (also known as DNp37, drives vaginal plate opening, a receptivity behavior) \- if any of these connected cell types are in MANC, then make a comparison \- 2 of these DN types are female-specific (oviDNs and vpoDN) and one is sex-dimorphic (DNp13, which is more involved in song production in males)) 

**Circuitry Downstream of DNp13 in Males and Females**  
	Previous work suggests that DNp13 has major outputs in the wing tectulum in the female (Sturner, IN06B035, IN06B047, IN06B050). While we do observe wing output overlapping with this report at the \>=2% levels (female only edges: IN17A114, IN06B050, IN06B038, IN06B071, IN17A119, IN23B012, IN08B104; edges shared with males: IN06B047, IN03B057, IN12A002), the major outputs of female DNp13 are concentrated in the AbG. 

FIgure 6 (if time): Sex Dimorphisms in Forward and Backward Walking:  
DNg100 (aka BDN2) drives forward walking and MDN drives backward walking. Compare how these circuits downstream of these two DN types differ in males and females

**Catalog of Sex Differences** (presented in Supplemental Materials?):

Here provide an overview of resources available in Codex for exploring the sex different circuitry in the VNCs of BANC and MANC \- and how people can use the catalog for future science

# **References**

[Bates, Alexander Shakeel, Jasper S. Phelps, Minsu Kim, et al. 2025\. “Distributed Control Circuits across a Brain-and-Cord Connectome.” In *bioRxivorg*. August 2\. https://doi.org/](http://paperpile.com/b/7p9NJj/ETyZ)[10.1101/2025.07.31.667571](http://dx.doi.org/10.1101/2025.07.31.667571)[.](http://paperpile.com/b/7p9NJj/ETyZ)

[Berg, Stuart, Isabella R. Beckett, Marta Costa, et al. 2025\. “Sexual Dimorphism in the Complete Connectome of the *Drosophila* Male Central Nervous System.” In *bioRxiv*. October 9\. https://doi.org/](http://paperpile.com/b/7p9NJj/aKLi)[10.1101/2025.10.09.680999](http://dx.doi.org/10.1101/2025.10.09.680999)[.](http://paperpile.com/b/7p9NJj/aKLi)

[Costa, Marta, James D. Manton, Aaron D. Ostrovsky, Steffen Prohaska, and Gregory S. X. E. Jefferis. 2016\. “NBLAST: Rapid, Sensitive Comparison of Neuronal Structure and Construction of Neuron Family Databases.” *Neuron* 91 (2): 293–311.](http://paperpile.com/b/7p9NJj/15Ad)

[Takemura, Shin-Ya, Kenneth J. Hayworth, Gary B. Huang, et al. 2023\. “A Connectome of the male*Drosophila*ventral Nerve Cord.” In *bioRxiv*. June 6\. https://doi.org/](http://paperpile.com/b/7p9NJj/E9X8)[10.1101/2023.06.05.543757](http://dx.doi.org/10.1101/2023.06.05.543757)[.](http://paperpile.com/b/7p9NJj/E9X8)

Useful review? [https://www.nature.com/articles/s41592-025-02946-2](https://www.nature.com/articles/s41592-025-02946-2)

1. Connectome comparison methods: [https://www.nature.com/articles/s41586-024-07686-5](https://www.nature.com/articles/s41586-024-07686-5); [https://www.biorxiv.org/content/10.1101/2025.07.16.664682v3.full.pdf](https://www.biorxiv.org/content/10.1101/2025.07.16.664682v3.full.pdf); [https://www.biorxiv.org/content/10.1101/2025.10.09.680999v2](https://www.biorxiv.org/content/10.1101/2025.10.09.680999v2)  
2. Drosophila sex dimorphism (dsx/fru): [https://www.biorxiv.org/content/10.1101/2025.06.04.657833v1](https://www.biorxiv.org/content/10.1101/2025.06.04.657833v1); [https://www.biorxiv.org/content/10.1101/2025.06.10.658788v1](https://www.biorxiv.org/content/10.1101/2025.06.10.658788v1); [https://www.sciencedirect.com/science/article/pii/S0960982210009474](https://www.sciencedirect.com/science/article/pii/S0960982210009474)  
3. Brain connectomics:  
4. VNC connectomics: MANC: [https://elifesciences.org/reviewed-preprints/97769](https://elifesciences.org/reviewed-preprints/97769); MANC: [https://elifesciences.org/reviewed-preprints/97766](https://elifesciences.org/reviewed-preprints/97766); MANC: [https://elifesciences.org/reviewed-preprints/96084](https://elifesciences.org/reviewed-preprints/96084); FANC: [https://www.nature.com/articles/s41586-024-07600-z](https://www.nature.com/articles/s41586-024-07600-z); FANC: [https://www.cell.com/cell/fulltext/S0092-8674(20)31683-4?uuid=uuid%3Ab0357d9b-f22c-4fc5-a203-63138d1d9b9d](https://www.cell.com/cell/fulltext/S0092-8674\(20\)31683-4?uuid=uuid%3Ab0357d9b-f22c-4fc5-a203-63138d1d9b9d)  
5. VNC sex dimorphisms and transcriptomes: [https://www.biorxiv.org/content/10.1101/2025.07.16.664682v3.full.pdf](https://www.biorxiv.org/content/10.1101/2025.07.16.664682v3.full.pdf); [https://elifesciences.org/articles/54074](https://elifesciences.org/articles/54074)  
6. Anatomy of the VNC: [https://www.sciencedirect.com/science/article/pii/S0896627320306127](https://www.sciencedirect.com/science/article/pii/S0896627320306127);   
7. Functional imaging of VNC: [https://www.nature.com/articles/s41467-022-32571-y](https://www.nature.com/articles/s41467-022-32571-y); [https://www.nature.com/articles/s41467-018-06857-z](https://www.nature.com/articles/s41467-018-06857-z)  
8. DN connectomics and cell types: [https://www.nature.com/articles/s41586-025-08925-z](https://www.nature.com/articles/s41586-025-08925-z); [https://elifesciences.org/articles/34272](https://elifesciences.org/articles/34272)  
9. Male song circuits (melanogaster): [https://www.sciencedirect.com/science/article/pii/S0960982224000150](https://www.sciencedirect.com/science/article/pii/S0960982224000150); [https://pubmed.ncbi.nlm.nih.gov/27326931/](https://pubmed.ncbi.nlm.nih.gov/27326931/); [https://www.nature.com/articles/s41593-024-01738-9](https://www.nature.com/articles/s41593-024-01738-9); [https://www.cell.com/cell-reports/pdfExtended/S2211-1247(13)00562-7](https://www.cell.com/cell-reports/pdfExtended/S2211-1247\(13\)00562-7); [https://www.sciencedirect.com/science/article/pii/S0896627311000572](https://www.sciencedirect.com/science/article/pii/S0896627311000572); [https://www.nature.com/articles/s41586-023-06632-1](https://www.nature.com/articles/s41586-023-06632-1); [https://www.sciencedirect.com/science/article/pii/S0960982218307735](https://www.sciencedirect.com/science/article/pii/S0960982218307735);   
10. Male song circuits (non-melanogaster): [https://www.cell.com/current-biology/fulltext/S0960-9822(24)00460-3?rss=yes](https://www.cell.com/current-biology/fulltext/S0960-9822\(24\)00460-3?rss=yes)**;**   
11. vpoDN (DNp37): [https://pubmed.ncbi.nlm.nih.gov/33239786/](https://pubmed.ncbi.nlm.nih.gov/33239786/); [https://www.nature.com/articles/s41467-024-53610-w](https://www.nature.com/articles/s41467-024-53610-w)  
12. DNp13: [https://www.cell.com/current-biology/pdfExtended/S0960-9822(20)31142-8](https://www.cell.com/current-biology/pdfExtended/S0960-9822\(20\)31142-8); [https://www.cell.com/current-biology/fulltext/S0960-9822(20)30923-4](https://www.cell.com/current-biology/fulltext/S0960-9822\(20\)30923-4)  
13. oviDN: [https://www.nature.com/articles/s41586-020-2055-9](https://www.nature.com/articles/s41586-020-2055-9); [https://www.nature.com/articles/s41586-023-06271-6](https://www.nature.com/articles/s41586-023-06271-6); [https://www.science.org/doi/10.1126/sciadv.abn3852](https://www.science.org/doi/10.1126/sciadv.abn3852); [https://www.nature.com/articles/s41593-023-01332-5](https://www.nature.com/articles/s41593-023-01332-5)  
14. Muscle and motor neuron activity during flight and song: [https://www.sciencedirect.com/science/article/pii/S0960982216314658](https://www.sciencedirect.com/science/article/pii/S0960982216314658); [https://www.sciencedirect.com/science/article/pii/S0960982218308297](https://www.sciencedirect.com/science/article/pii/S0960982218308297); [https://elifesciences.org/reviewed-preprints/106548](https://elifesciences.org/reviewed-preprints/106548); [https://journals.biologists.com/jeb/article-abstract/49/1/117/21364/The-Wing-Mechanism-Involved-in-the-Courtship-of?redirectedFrom=fulltext](https://journals.biologists.com/jeb/article-abstract/49/1/117/21364/The-Wing-Mechanism-Involved-in-the-Courtship-of?redirectedFrom=fulltext); [https://resjournals.onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-3032.1979.tb00624.x](https://resjournals.onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-3032.1979.tb00624.x)

Figures 1 and 2:

1. Method   
2. Test on maleCNS VNC to MANC matching \- establishing baseline for sex differences (versus individual differences)

TODO for data:

1. Finish typing of BANC Sex Specific and Dimorphic ANs  
2. Finish assembling untyped cells into types  
3. Rectify types  
   1. Remove dissimilar cells from types and reassign  
4. Find matches for unmatched homologous types (see mismatches in banc\_manc\_descriptive.R)  
5. Warp BANC and MANC into JRC2018 space for synapse density and cell morpholohy comparison

Chris TODO for Figs:

1. Descriptive quantification of datasets by Neuromere:   
   1. cell numbers  
   2. cell types  
   3. cells per type  
   4. Sex specific and dimorphic cells. Secondarily dimorphic as well?  
   5. **I am porting type neuromere assignments to BANC from MANC for this**  
2. Distribution of BANC-MANC cell type similarity  
   1. This should be able to come right from the inventory page I think?  
   2. Overlay the distribution of dimorphic and sex specific cells on the same graph  
   3. Also do it separately for  
      1. IN  
      2. MN  
      3. AN  
      4. DN  
      5. SN  
3. Examples of new dimorphic types here?  
   1. Wing  
   2. Legs  
   3. AbG  
4. Overview Image of synapses made by sex-specific and dimorphic neurons  
   1. Both BANC and MANC warped into JRC2018 space  
   2. Heat map for synapse density  
5. Quantification of synapses made by sex specific and dimorphic cells in different Neuromeres and ideally Neuropils if we can get that info  
   1. For dimorphs, then break it down by sex-specific and homologous edges  
      1. **Arie/Alex \- This will require neuromere and neuropil columns in BANC [Arie Matsliah](mailto:am1736@princeton.edu)**  
   2. For sex specific, break it down by synapses made onto sex specific, dimorphic and homologous neurons.  
      1. **Arie/Alex \- this will require neuromere and neuropil columns for BANC and MANC *edges*** [Arie Matsliah](mailto:am1736@princeton.edu)  
6. 
