# NMF analysis


## Applying NMF
After preprocessing the data, we are now ready to start an NMF analysis.

------------

### Call wrapper function

The wrapper function for the NMF solvers in the ButchR package is
`run_NMF_tensor`. It is called as follows:


Applying Non-Negative Matrix Factorization (NMF) to normalized transcriptome data (RNA-seq). Using:  
- Factorization ranks from 2 to 10.
- Default NMF method.
- 10 random initialization.

Depending on the choice of parameters (dimensions of the input matrix,
number of iterations), this step may take some time. Note that the
algorithm updates the user about the progress in the iterations.

```r
##----------------------------------------------------------------------------##
##                             run NMF                                        ##
##----------------------------------------------------------------------------##
factorization_ranks <- 2:10
rna_nmf_exp <- run_NMF_tensor(X                     = corces_rna_norm,
                              ranks                 = factorization_ranks,
                              method                = "NMF",
                              n_initializations     = 10,
                              iterations            = 10^4,
                              convergence_threshold = 40, 
                              extract_features = TRUE)
```

<details>
<summary><b>Click for Answer</b></summary>
```
## [1] "2021-09-02 16:23:14 UTC"
## Factorization rank:  2 
## [1] "NMF converged after  71,106,122,177,111,84,168,66,84,106 iterations"
## [1] "2021-09-02 16:23:19 UTC"
## Factorization rank:  3 
## [1] "NMF converged after  141,152,77,129,177,112,95,159,118,91 iterations"
## [1] "2021-09-02 16:23:27 UTC"
## Factorization rank:  4 
## [1] "NMF converged after  122,244,107,213,94,138,239,138,140,93 iterations"
## [1] "2021-09-02 16:23:36 UTC"
## Factorization rank:  5 
## [1] "NMF converged after  182,187,151,174,310,244,192,164,187,191 iterations"
## [1] "2021-09-02 16:23:47 UTC"
## Factorization rank:  6 
## [1] "NMF converged after  159,148,142,215,208,166,153,273,267,355 iterations"
## [1] "2021-09-02 16:23:59 UTC"
## Factorization rank:  7 
## [1] "NMF converged after  220,171,205,313,116,249,168,163,145,295 iterations"
## [1] "2021-09-02 16:24:13 UTC"
## Factorization rank:  8 
## [1] "NMF converged after  217,176,199,289,257,127,216,153,132,138 iterations"
## [1] "2021-09-02 16:24:26 UTC"
## Factorization rank:  9 
## [1] "NMF converged after  167,209,176,155,186,253,241,136,161,191 iterations"
## [1] "2021-09-02 16:24:39 UTC"
## Factorization rank:  10 
## [1] "NMF converged after  304,265,312,368,301,214,172,295,161,279 iterations"
## No optimal K could be determined from the Optimal K stat
```
</details>

```r
rna_nmf_exp
```

<details>
<summary><b>Click for Answer</b></summary>
```
## class: ButchR_NMF 
## Original matrix dimension:  21811 45 
## Factorization performed for ranks:  2 3 4 5 6 7 8 9 10 
## Optimal K based on factorization metrics:  Please select manualy
##  
## Running parameters: 
## method =  NMF  
## n_initializations =  10  
## iterations =  10000  
## stop_threshold =  40  
## extract_features =  TRUE
```
</details>



### Normalize W matrix

To make the features in the *W* matrix comparable, the factorization is
normalized to make all columns of *W* sum 1.

```r
rna_norm_nmf_exp <- normalizeW(rna_nmf_exp)
```
 
## Accessor functions

Several functions to access the results are available:

### `HMatrix`

Returns the matrix `H` for the optimal decomposition (i.e. the one with
the minimal residual) for a specific factorization rank `k`. The number
of rows of the matrix `H` corresponds to the chosen factorization rank.

Extract the matrix H for k=2 and find its dimension:
```r
Hk2 <- HMatrix(rna_norm_nmf_exp, k = 2)
class(Hk2)
```
<details>
<summary><b>Click for Answer</b></summary>
```
## [1] "matrix" "array"
```

```r
dim(Hk2)
```

```
## [1]  2 45
```
</details>

Inspect exposure values of the first 5 samples in the H matrix k=2:

```r
Hk2[, 1:5]
```
<details>
<summary><b>Click for Answer</b></summary>
<table>
 <thead>
  <tr>
   <th style="text-align:right;"> X5852.HSC </th>
   <th style="text-align:right;"> X6792.HSC </th>
   <th style="text-align:right;"> X7256.HSC </th>
   <th style="text-align:right;"> X7653.HSC </th>
   <th style="text-align:right;"> X5852.MPP </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 73737.72 </td>
   <td style="text-align:right;"> 73937.34 </td>
   <td style="text-align:right;"> 66949.25 </td>
   <td style="text-align:right;"> 74760.02 </td>
   <td style="text-align:right;"> 76202.74 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 20470.68 </td>
   <td style="text-align:right;"> 21970.29 </td>
   <td style="text-align:right;"> 26204.85 </td>
   <td style="text-align:right;"> 20569.97 </td>
   <td style="text-align:right;"> 18053.79 </td>
  </tr>
</tbody>
</table>
</details>


If no value for `k` is supplied, the function returns a list of
matrices, one for every factorization rank.

```r
Hlist <- HMatrix(rna_norm_nmf_exp)
class(Hlist)
```

<details>
<summary><b>Click for Answer</b></summary>
```
## [1] "list"
```

```r
length(Hlist)
```

```
## [1] 9
```
</details>


```r
Hlist$k2[, 1:5]
```

<details>
<summary><b>Click for Answer</b></summary>
<table>
 <thead>
  <tr>
   <th style="text-align:right;"> X5852.HSC </th>
   <th style="text-align:right;"> X6792.HSC </th>
   <th style="text-align:right;"> X7256.HSC </th>
   <th style="text-align:right;"> X7653.HSC </th>
   <th style="text-align:right;"> X5852.MPP </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 73737.72 </td>
   <td style="text-align:right;"> 73937.34 </td>
   <td style="text-align:right;"> 66949.25 </td>
   <td style="text-align:right;"> 74760.02 </td>
   <td style="text-align:right;"> 76202.74 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 20470.68 </td>
   <td style="text-align:right;"> 21970.29 </td>
   <td style="text-align:right;"> 26204.85 </td>
   <td style="text-align:right;"> 20569.97 </td>
   <td style="text-align:right;"> 18053.79 </td>
  </tr>
</tbody>
</table>
</details>


### `WMatrix`

Returns the matrix `W` for the optimal decomposition (i.e. the one with
the minimal residual) for a specific factorization rank `k`. The number
of columns of the matrix `W` corresponds to the chosen factorization
rank.

Extract the matrix W for k=2 and find its dimension:

```r
Wk2 <- WMatrix(rna_norm_nmf_exp, k = 2)
class(Wk2)
```

<details>
<summary><b>Click for Answer</b></summary>
```
## [1] "matrix" "array"
```

    


```r
dim(Wk2)
```

```
## [1] 21811     2
```
</details>


Inspect the contribution of the first 5 genes to the signatures in matrix W k=2:
```r
Wk2[1:5, ]
```

<details>
<summary><b>Click for Answer</b></summary>
<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> V1 </th>
   <th style="text-align:right;"> V2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> A1BG </td>
   <td style="text-align:right;"> 1.99e-05 </td>
   <td style="text-align:right;"> 2.00e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A1BG-AS1 </td>
   <td style="text-align:right;"> 1.21e-05 </td>
   <td style="text-align:right;"> 9.50e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A1CF </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 1.34e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A2M </td>
   <td style="text-align:right;"> 3.88e-05 </td>
   <td style="text-align:right;"> 8.75e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A2M-AS1 </td>
   <td style="text-align:right;"> 5.25e-05 </td>
   <td style="text-align:right;"> 4.86e-05 </td>
  </tr>
</tbody>
</table>
</details>


If no value for `k` is supplied, the function returns a list of
matrices, one for every factorization rank.


```r
Wlist <- WMatrix(rna_norm_nmf_exp)
class(Wlist)
```

<details>
<summary><b>Click for Answer</b></summary>
```
## [1] "list"
```

    

```r
length(Wlist)
```

```
## [1] 9
```
</details>

    
```r
Wlist$k2[1:5, ]
```

<details>
<summary><b>Click for Answer</b></summary>
<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> V1 </th>
   <th style="text-align:right;"> V2 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> A1BG </td>
   <td style="text-align:right;"> 1.99e-05 </td>
   <td style="text-align:right;"> 2.00e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A1BG-AS1 </td>
   <td style="text-align:right;"> 1.21e-05 </td>
   <td style="text-align:right;"> 9.50e-06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A1CF </td>
   <td style="text-align:right;"> 0.00e+00 </td>
   <td style="text-align:right;"> 1.34e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A2M </td>
   <td style="text-align:right;"> 3.88e-05 </td>
   <td style="text-align:right;"> 8.75e-05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> A2M-AS1 </td>
   <td style="text-align:right;"> 5.25e-05 </td>
   <td style="text-align:right;"> 4.86e-05 </td>
  </tr>
</tbody>
</table>
</details>
