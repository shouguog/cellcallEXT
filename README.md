<a name="PNqrO"></a>
# CellCallEXT: An extention of integrative analysis of paired ligand–receptor and transcription factor activities for cell–cell communication 


<a name="adf749d9"></a>
## Updated information of CellCall

<a name="c8f560dd"></a>
## 1. Introduction to CellCallEXT

<a name="JqFcu"></a>
### 1.1 Overview of CellCallEXT
CellCallEXT extended R package for cellcall, with large databases and enabling to analysis case control dataset.
### 1.2 Installing R package
To install this package, start R (version 3.6 or higher) and enter:
```
library(devtools)
devtools::install_github("ShellyCoder/cellcall")
```
If you encounter the following error -- ERROR: dependency is not available for package 'cellcall', install corresponding R package. And appropriate version is in the section 4.<br />

<a name="689d6784"></a>
## 2. Main functions of CellCall
CellCall provides a variety of functions including intercellular communication analysis, pathway activity analysis and a rich suite of visualization tools to intuitively present the results of the analysis (including Heatmap, Circos plot, Bubble plot, Sankey plot, TF enrichment plot and Ridge plot).<br />

<a name="LP04w"></a>
### 2.1 Intercellular communication analysis


<a name="945a6b52"></a>
#### 2.1.1 Load data
The format of the input file is as follow table:<br />1. The row names: gene symbols.<br />2. The column names: cell IDs. The colnames can't contain punctuation such as commas, periods, dashes, etc. Using underline to connect barcoder and cell type is recommended. Take the input format below as an example, the column name is made up of index and cell type. Users should set names.field=2 and names.delim="_" in the function CreateNichConObject(). After that, cell type information is obtained and stored in the S4 object for later analysis. Because method in this paper depends on the cell type information, obtaining celltype information correctly is important.<br />3. Other place: the expression values (counts or TPM) for a gene in a cell.

|  | 1_ST | 2_ST | 3_ST | 4_SSC | 5_SSC | 6_SPGed | 7_SPGed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| TSPAN6 | 2.278 | 2.031 | 0.000 | 12.385 | 0.000 | 0.553 | 24.846 |
| TNMD | 9.112 | 6.031 | 0.000 | 0.000 | 11.615 | 10.518 | 0.000 |
| DPM1 | 0.000 | 0.000 | 21.498 | 4.246 | 7.382 | 0.000 | 2.385 |
| SCYL3 | 5.983 | 1.215 | 0.000 | 0.518 | 2.386 | 4.002 | 14.792 |

This instruction may take the in-house dataset included in the package as an example. User can load the dataset with command following in the code box. There are 366 single cells and 35,135 genes that were performed with the scRNA sequencing.
```r
f.tmp <- system.file("extdata", "example_Data.Rdata", package="cellcallEXT")
load(f.tmp)

## gene expression stored in the variable in.content
dim(in.content)
in.content[1:4, 1:4]
table(str_split(colnames(in.content), "_", simplify = T)[,2])
```
We next use the expression dataframe  to create a CreateNichCon object with the function CreateNichConObject as the code in part 2.1.2. The object serves as a container that contains both data (like the expression dataframe) and analysis (like score, or enrichment results) for a single-cell dataset.<br />

<a name="e815f697"></a>
#### 2.1.2 Create object
CellCall use the expression dataframe to create an S4 object by the function CreateNichConObject, The line of code is shown in the code box. The S4 object serves as a container that contains both data (such as the expression dataframe) and analysis results (such score and enrichment results) for a project (see section 3 for details).
```r
in.content[1:4, 1:4]; set.seed(666); in.content<-in.content[, sample(1:366, 100)]
set.seed(666)
mt <- CreateNichConObject(data=in.content, min.feature = 3,
                            names.field = 2,
                            names.delim = "_",
                            source = "TPM",
                            scale.factor = 10^6,
                            Org = "Homo sapiens",
                            status = c(rep("CASE", 50), rep("CONTROL", 50))[sample(1:100, 100)],
                            expDirection = "UP"
                            project = "Microenvironment")
```
**Arguments**:

| **Arguments** | **Detail** |
| --- | --- |
| **data** | A dataframe with row of gene and column of sample and the value must be numeric. Meanwhile what matters is that the colnames of dataframe should be in line with the paramter 'names.delim' and 'names.field', the former for pattern to splite every colnames, the latter for setting which index in splited colnames is cell type information.<br />The function can get the 'CELLTYPE' information from the colnames 'BARCODE_CLUSTER_CELLTYPE' with names.delim="_" and names.field='3', and then stored in slot meta.data of CreateNichCon.<br />Cell type annotation from every cell is essential for scoring cell communication. If the colnames of data don't coincide with the paramter 'names.delim' and 'names.field', CreateNichCon object may fail to create. |
| **min.feature** | Include cells where enough features equalling min.feature are detected. It's a preprocess which is the same as Seurat and set min.feature=0, if you don't want to filter cell. This parameter depends on the sequencing technology of the input data. |
| **names.delim** | Set the pattern to splite column name into vector. If the column name of the input matrix is BARCODE_CLUSTER_CELLTYPE, set names.delim="_" to get CELLTYPE of BARCODE_CLUSTER_CELLTYPE with names.field=3. |
| **names.field** | Set the index of column name vector which is splited by parameter names.delim to get cell type information. If the column name of the input matrix is BARCODE_CLUSTER_CELLTYPE, set names.field=3 to get CELLTYPE of BARCODE_CLUSTER_CELLTYPE with names.delim="_". |
| **source** | The type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM". When the source of input data is  "TPM" or "CPM", no transformation on the data. Otherwise, we transform the data to TPM with the parameter source="fullLength" and to CPM with source="UMI". |
| **scale.factor** | Sets the scale factor for cell-level normalization, default "10^6", if the parameter is "UMI" or "fullLength". Otherwise this parameter doesn't work. |
| **Org** | Set the species source of gene, eg "Homo sapiens", "Mus musculus". This parameter decides the paired ligand-receptor dataset and the transcript length which is needed in "TPM" transformation. |
| **status** | the status of cell disease or control. |
| **expDirection** | the direction of expression change UP DOWN or BOTH. |
| **project** | Sets the project name for the CreateNichCon object. |


<br />

<a name="oEar5"></a>
#### 2.1.3 Infer the cell-cell communication score
The communication score of an L-R interaction between cell types is evaluated by integrating the L2- norm of the L-R interaction and the activity score of the downstream TFs. The code is shown in the code box.
```r
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Homo sapiens')
```
