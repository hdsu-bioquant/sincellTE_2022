---
title: "IRTG Course - Introduction to R for Genomics"
author: "Carl Herrmann & Carlos Ramirez"
date: "2021-06-12"
output: 
  html_document: 
    keep_md: yes
---





# 3. Data filtering and cleanup

Very often the first thing one needs to do before any data science project is to clean up the raw data and transform it into a format that is readily understood and easy to use for all downstream analysis. This process usually involves:

* Removing empty value rows/columns
* Removing unused or unnecessary rows/columns
* Reordering the data matrix
* Keeping columns uniformly numeric (age, weight etc) or string (names, places etc) or logical (TRUE/FALSE, 1/0)
* Handling strange caveats which are data specific like replacing `,` or `.`, or `;` from numbers etc

In addition, you often want to filter rows according to vertain criteria: for example, selecting only women of age more than 45.

Here, we will learn how to

* clean up the data
* filter the data

In order to do these manipulations easily, we will rely on a library which has a lot of functions to easily manipulate tables: [dplyr](https://dplyr.tidyverse.org/). This library is part of a large eco-system of data analysis called **tidyverse**.

Let's load the library


```r
library(tidyverse)
```

We read again our dataset:


```r
dat = read.delim("https://tinyurl.com/ex9hxvvr", stringsAsFactors = FALSE)
```


Lets do some clean up of our own diabetes data --

1. We will make the `id` column the row names for the dataset. 
2. We will remove the `bp.2s` and `bp.2d` columns as it has mostly missing values (remember the output of the `summary` function in part 1? If not, run it again...)
3. We will also remove the column `time.ppn` which will not be required in our analysis
4. We will reorder the columns of the data such that all the qualitative and quantitative values are separated. Among the quantitative values we will keep related variables together

## Setting row names

As it is, our data frame has no row names (check by running `rownames(dat)`); however, it might be interesting to have row names.

**IMPORTANT**: for data frames, all row names must be distinct!! This can be a problem if you are working with genes, which have sometimes ambiguous names. 
Here, we will assign the column `id` as row names:


```r
## Assign the id column as row names
dat.clean = dat %>% column_to_rownames(c("id"))
```

<details>
<summary><b>More information: the pipe function %>%</b></summary>

You noticed how we used the syntax `dat %>% column_to_rownames(c('id'))`; what does it mean?

The symbol `%>%` is a so called pipe, and as the name indicates, it "pipes" the content on the left hand side to the function on the right hand side.

Here, the content of the variable `dat` is transfered to the function `column_to_rownames(c('id'))` which then applies the operation (here, using column `id` as rownames, and then removing this column) to the content of the data frame.

Note that we also assigned the result of this whole operation to a new variable `data.clean`.

We could extend the pipe principle and chain several commands!
</details>
<p></p>

## Reordering columns

Now, we want to reorder the columns so that all numerical columns come first, and the other columns after:

The `relocate` function allows to put a column to a certain position: try the following:


```r
head(dat.clean)
```

```
     chol stab.glu hdl ratio glyhb   location age gender height weight  frame
1000  203       82  56   3.6  4.31 Buckingham  46 female     62    121 medium
1001  165       97  24   6.9  4.44 Buckingham  29 female     64    218  large
1002  228       92  37   6.2  4.64 Buckingham  58 female     61    256  large
1003   78       93  12   6.5  4.63 Buckingham  67   male     67    119  large
1005  249       90  28   8.9  7.72 Buckingham  64   male     68    183 medium
1008  248       94  69   3.6  4.81 Buckingham  34   male     71    190  large
     bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1000   118    59    NA    NA    29  38      720
1001   112    68    NA    NA    46  48      360
1002   190    92   185    92    49  57      180
1003   110    50    NA    NA    33  38      480
1005   138    80    NA    NA    44  41      300
1008   132    86    NA    NA    36  42      195
```

Notice the order of the columns!
Now:


```r
dat.clean %>% relocate(gender)
```

```
     gender chol stab.glu hdl ratio glyhb   location age height weight  frame
1000 female  203       82  56   3.6  4.31 Buckingham  46     62    121 medium
1001 female  165       97  24   6.9  4.44 Buckingham  29     64    218  large
1002 female  228       92  37   6.2  4.64 Buckingham  58     61    256  large
1003   male   78       93  12   6.5  4.63 Buckingham  67     67    119  large
1005   male  249       90  28   8.9  7.72 Buckingham  64     68    183 medium
1008   male  248       94  69   3.6  4.81 Buckingham  34     71    190  large
1011   male  195       92  41   4.8  4.84 Buckingham  30     69    191 medium
1015   male  227       75  44   5.2  3.94 Buckingham  37     59    170 medium
1016   male  177       87  49   3.6  4.84 Buckingham  45     69    166  large
1022 female  263       89  40   6.6  5.78 Buckingham  55     63    202  small
1024 female  242       82  54   4.5  4.77     Louisa  60     65    156 medium
     bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1000   118    59    NA    NA    29  38      720
1001   112    68    NA    NA    46  48      360
1002   190    92   185    92    49  57      180
1003   110    50    NA    NA    33  38      480
1005   138    80    NA    NA    44  41      300
1008   132    86    NA    NA    36  42      195
1011   161   112   161   112    46  49      720
1015    NA    NA    NA    NA    34  39     1020
1016   160    80   128    86    34  40      300
1022   108    72    NA    NA    45  50      240
1024   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

Notice the difference? You can also position a column to a specific position, before or after a specified column:


```r
dat.clean %>% relocate(gender, .before = age)
```

```
     chol stab.glu hdl ratio glyhb   location gender age height weight  frame
1000  203       82  56   3.6  4.31 Buckingham female  46     62    121 medium
1001  165       97  24   6.9  4.44 Buckingham female  29     64    218  large
1002  228       92  37   6.2  4.64 Buckingham female  58     61    256  large
1003   78       93  12   6.5  4.63 Buckingham   male  67     67    119  large
1005  249       90  28   8.9  7.72 Buckingham   male  64     68    183 medium
1008  248       94  69   3.6  4.81 Buckingham   male  34     71    190  large
1011  195       92  41   4.8  4.84 Buckingham   male  30     69    191 medium
1015  227       75  44   5.2  3.94 Buckingham   male  37     59    170 medium
1016  177       87  49   3.6  4.84 Buckingham   male  45     69    166  large
1022  263       89  40   6.6  5.78 Buckingham female  55     63    202  small
1024  242       82  54   4.5  4.77     Louisa female  60     65    156 medium
     bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1000   118    59    NA    NA    29  38      720
1001   112    68    NA    NA    46  48      360
1002   190    92   185    92    49  57      180
1003   110    50    NA    NA    33  38      480
1005   138    80    NA    NA    44  41      300
1008   132    86    NA    NA    36  42      195
1011   161   112   161   112    46  49      720
1015    NA    NA    NA    NA    34  39     1020
1016   160    80   128    86    34  40      300
1022   108    72    NA    NA    45  50      240
1024   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

Here, we would like to place all numeric columns first, and the columns containing string or other types of variables after:


```r
dat.clean = dat.clean %>% relocate(where(is.numeric))
head(dat.clean)
```

```
     chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d bp.2s bp.2d
1000  203       82  56   3.6  4.31  46     62    121   118    59    NA    NA
1001  165       97  24   6.9  4.44  29     64    218   112    68    NA    NA
1002  228       92  37   6.2  4.64  58     61    256   190    92   185    92
1003   78       93  12   6.5  4.63  67     67    119   110    50    NA    NA
1005  249       90  28   8.9  7.72  64     68    183   138    80    NA    NA
1008  248       94  69   3.6  4.81  34     71    190   132    86    NA    NA
     waist hip time.ppn   location gender  frame
1000    29  38      720 Buckingham female medium
1001    46  48      360 Buckingham female  large
1002    49  57      180 Buckingham female  large
1003    33  38      480 Buckingham   male  large
1005    44  41      300 Buckingham   male medium
1008    36  42      195 Buckingham   male  large
```


> try to place all string variables first; use the `is.character` function!

<details>
<summary><b>Click for solution!</b></summary> 


```r
dat.clean %>% relocate(where(is.character))
```

```
       location gender  frame chol stab.glu hdl ratio glyhb age height weight
1000 Buckingham female medium  203       82  56   3.6  4.31  46     62    121
1001 Buckingham female  large  165       97  24   6.9  4.44  29     64    218
1002 Buckingham female  large  228       92  37   6.2  4.64  58     61    256
1003 Buckingham   male  large   78       93  12   6.5  4.63  67     67    119
1005 Buckingham   male medium  249       90  28   8.9  7.72  64     68    183
1008 Buckingham   male  large  248       94  69   3.6  4.81  34     71    190
1011 Buckingham   male medium  195       92  41   4.8  4.84  30     69    191
1015 Buckingham   male medium  227       75  44   5.2  3.94  37     59    170
1016 Buckingham   male  large  177       87  49   3.6  4.84  45     69    166
1022 Buckingham female  small  263       89  40   6.6  5.78  55     63    202
1024     Louisa female medium  242       82  54   4.5  4.77  60     65    156
     bp.1s bp.1d bp.2s bp.2d waist hip time.ppn
1000   118    59    NA    NA    29  38      720
1001   112    68    NA    NA    46  48      360
1002   190    92   185    92    49  57      180
1003   110    50    NA    NA    33  38      480
1005   138    80    NA    NA    44  41      300
1008   132    86    NA    NA    36  42      195
1011   161   112   161   112    46  49      720
1015    NA    NA    NA    NA    34  39     1020
1016   160    80   128    86    34  40      300
1022   108    72    NA    NA    45  50      240
1024   130    90   130    90    39  45      300
 [ reached 'max' / getOption("max.print") -- omitted 392 rows ]
```

</details> 
<p></p>

Now lets look at our cleaned data

```r
summary(dat.clean)
```

```
      chol          stab.glu          hdl             ratio       
 Min.   : 78.0   Min.   : 48.0   Min.   : 12.00   Min.   : 1.500  
 1st Qu.:179.0   1st Qu.: 81.0   1st Qu.: 38.00   1st Qu.: 3.200  
 Median :204.0   Median : 89.0   Median : 46.00   Median : 4.200  
 Mean   :207.8   Mean   :106.7   Mean   : 50.45   Mean   : 4.522  
 3rd Qu.:230.0   3rd Qu.:106.0   3rd Qu.: 59.00   3rd Qu.: 5.400  
 Max.   :443.0   Max.   :385.0   Max.   :120.00   Max.   :19.300  
 NA's   :1                       NA's   :1        NA's   :1       
     glyhb            age            height          weight     
 Min.   : 2.68   Min.   :19.00   Min.   :52.00   Min.   : 99.0  
 1st Qu.: 4.38   1st Qu.:34.00   1st Qu.:63.00   1st Qu.:151.0  
 Median : 4.84   Median :45.00   Median :66.00   Median :172.5  
 Mean   : 5.59   Mean   :46.85   Mean   :66.02   Mean   :177.6  
 3rd Qu.: 5.60   3rd Qu.:60.00   3rd Qu.:69.00   3rd Qu.:200.0  
 Max.   :16.11   Max.   :92.00   Max.   :76.00   Max.   :325.0  
 NA's   :13                      NA's   :5       NA's   :1      
     bp.1s           bp.1d            bp.2s           bp.2d       
 Min.   : 90.0   Min.   : 48.00   Min.   :110.0   Min.   : 60.00  
 1st Qu.:121.2   1st Qu.: 75.00   1st Qu.:138.0   1st Qu.: 84.00  
 Median :136.0   Median : 82.00   Median :149.0   Median : 92.00  
 Mean   :136.9   Mean   : 83.32   Mean   :152.4   Mean   : 92.52  
 3rd Qu.:146.8   3rd Qu.: 90.00   3rd Qu.:161.0   3rd Qu.:100.00  
 Max.   :250.0   Max.   :124.00   Max.   :238.0   Max.   :124.00  
 NA's   :5       NA's   :5        NA's   :262     NA's   :262     
     waist           hip           time.ppn        location        
 Min.   :26.0   Min.   :30.00   Min.   :   5.0   Length:403        
 1st Qu.:33.0   1st Qu.:39.00   1st Qu.:  90.0   Class :character  
 Median :37.0   Median :42.00   Median : 240.0   Mode  :character  
 Mean   :37.9   Mean   :43.04   Mean   : 341.2                     
 3rd Qu.:41.0   3rd Qu.:46.00   3rd Qu.: 517.5                     
 Max.   :56.0   Max.   :64.00   Max.   :1560.0                     
 NA's   :2      NA's   :2       NA's   :3                          
    gender             frame          
 Length:403         Length:403        
 Class :character   Class :character  
 Mode  :character   Mode  :character  
                                      
                                      
                                      
                                      
```

The ordering and selection of columns looks right, however it seems that there are certain rows that have missing values (note which columns seem problematic!). 
We will now deal with the missing values.

## Dealing with NAs

Some columns and rows contain missing values, which are encoded by `NA` in the data frame. Missing values, especially if there are many can be an issue, as it will bias the results.
Dealing with missing values is a chapter in itself! Basically, there can be two strategies:

1. **imputing** missing value: this means that we try to make a "best guess" of what the missing value might be. For example, if the weight is missing for one patient, you could replace it with the average weight of the other patients. Of course, there are more sophisticated methods....

2. **removing** missing values: you could remove all patients (i.e. rows) are have
  + any missing value
  + more than a certain number of missing values
  + only missing values.
But you could also remove variables (i.e. columns) with a  lot of missing values!
  
Let's implement some of these strategies:

Before we do so, let's have a look at the `is.na` function


### Strategy 1: removing all rows which have **any** missing values:


```r
dat.nona = dat.clean %>% drop_na()
dim(dat.clean)
```

```
[1] 403  18
```

```r
dim(dat.nona)
```

```
[1] 130  18
```

> Do you understand why we reduced so dramatically our dataset?
> Could there be an alternative approach?

### Strategy 2: let us first remove problematic columns

Two of the columns have very high number of missing values `bp.2s` and `bp.2d`; we should remove these columns first and then remove rows with missing values:



```r
dat.clean = dat.clean %>% select(-c("bp.2s", "bp.2d"))
head(dat.clean)
```

```
     chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1000  203       82  56   3.6  4.31  46     62    121   118    59    29  38
1001  165       97  24   6.9  4.44  29     64    218   112    68    46  48
1002  228       92  37   6.2  4.64  58     61    256   190    92    49  57
1003   78       93  12   6.5  4.63  67     67    119   110    50    33  38
1005  249       90  28   8.9  7.72  64     68    183   138    80    44  41
1008  248       94  69   3.6  4.81  34     71    190   132    86    36  42
     time.ppn   location gender  frame
1000      720 Buckingham female medium
1001      360 Buckingham female  large
1002      180 Buckingham female  large
1003      480 Buckingham   male  large
1005      300 Buckingham   male medium
1008      195 Buckingham   male  large
```

We can now remove the rows with missing values:


```r
dat.nona = dat.clean %>% drop_na()
dim(dat.nona)
```

```
[1] 366  16
```

See? We have lost much less rows (or patients) here...

### strategy 3: remove rows which have more than a certain number of missing values

If you have many columns, it might be critical to remove rows with missing values in any column. Image a gene expression matrix with thousands of samples in columns: the likelihood that one gene might have a missing value in at least one of the samples is high!

Let's say that we want to remove rows which have more than 2 missing values (i.e. 3 or more); let's detail the procedure:

  + step 1: we need to determine the number of missing values in each row.


```r
dat.isna = is.na(dat.clean)
head(dat.isna)
```

```
      chol stab.glu   hdl ratio glyhb   age height weight bp.1s bp.1d waist
1000 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
1001 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
1002 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
1003 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
1005 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
1008 FALSE    FALSE FALSE FALSE FALSE FALSE  FALSE  FALSE FALSE FALSE FALSE
       hip time.ppn location gender frame
1000 FALSE    FALSE    FALSE  FALSE FALSE
1001 FALSE    FALSE    FALSE  FALSE FALSE
1002 FALSE    FALSE    FALSE  FALSE FALSE
1003 FALSE    FALSE    FALSE  FALSE FALSE
1005 FALSE    FALSE    FALSE  FALSE FALSE
1008 FALSE    FALSE    FALSE  FALSE FALSE
```

> What just happened? Do you understand how the function `is.na` works?

  + step 2: We now count the number of missing values in each row; we will sum up the number of `TRUE` in each row:



```r
nbr.na = rowSums(dat.isna)
head(nbr.na)
```

```
1000 1001 1002 1003 1005 1008 
   0    0    0    0    0    0 
```

  + step 3: We now use the `filter` function, which allows to filter rows according to certain criteria:
  

```r
dat.nona = dat.clean %>% filter(nbr.na < 3)
```
  
> How big is the `dat.nona` matrix?
> How many patients would you have kept if you had accept at most one missing value?

<details>
<summary><b>Click for solution!</b></summary> 


```r
dim(dat.nona)
```

```
[1] 400  16
```


```r
dim(dat.clean %>% filter(nbr.na < 2))
```

```
[1] 393  16
```
</details> 

Let's continue with the cleaned matrix obtained with strategy 2:



```r
dat.nona = dat.clean %>% drop_na()
```

## Converting strings into factors

Factors is a data type besides `numeric`, `characters` (or strings) `boolean`. It is very similar to the string type, but introduces the notion of **levels**, which indicates which categories are represented. Let's see an example:



```r
head(dat.nona$location)
```

```
[1] "Buckingham" "Buckingham" "Buckingham" "Buckingham" "Buckingham"
[6] "Buckingham"
```

These are strings. We cannot see easily however, how many different locations are represented in the column `location`. Let's convert this column into factors:


```r
dat.nona$location = factor(dat.nona$location)
head(dat.nona$location)
```

```
[1] Buckingham Buckingham Buckingham Buckingham Buckingham Buckingham
Levels: Buckingham Louisa
```

See the difference?

> Convert the `gender` and `frame` columns into factors!


<details>
<summary><b>Click for solution!</b></summary> 

```r
dat.nona$gender = factor(dat.nona$gender)  # Making data nominal
dat.nona$frame = factor(dat.nona$frame, levels = c("small", "medium", "large"))  # Making data ordinal
```
Notice that in the last command, we also indicated in which order the levels should be considered. Since we have ordinal data here (there is a clear order between small/medium/large), we should indicate this order here. Otherwise, the levels are ordered by alphabetical order!

</details> 
<p></p>

Let's inspect our cleaned dataset again:


```r
summary(dat.nona)
```

```
      chol          stab.glu          hdl             ratio       
 Min.   : 78.0   Min.   : 48.0   Min.   : 12.00   Min.   : 1.500  
 1st Qu.:179.0   1st Qu.: 81.0   1st Qu.: 38.00   1st Qu.: 3.200  
 Median :203.5   Median : 90.0   Median : 46.00   Median : 4.200  
 Mean   :207.5   Mean   :107.4   Mean   : 50.27   Mean   : 4.536  
 3rd Qu.:228.8   3rd Qu.:108.0   3rd Qu.: 59.00   3rd Qu.: 5.400  
 Max.   :443.0   Max.   :385.0   Max.   :120.00   Max.   :19.300  
     glyhb             age            height          weight     
 Min.   : 2.680   Min.   :19.00   Min.   :52.00   Min.   : 99.0  
 1st Qu.: 4.393   1st Qu.:34.00   1st Qu.:63.00   1st Qu.:151.0  
 Median : 4.860   Median :45.00   Median :66.00   Median :174.0  
 Mean   : 5.607   Mean   :46.69   Mean   :66.05   Mean   :178.1  
 3rd Qu.: 5.630   3rd Qu.:60.00   3rd Qu.:69.00   3rd Qu.:200.0  
 Max.   :16.110   Max.   :92.00   Max.   :76.00   Max.   :325.0  
     bp.1s           bp.1d            waist            hip       
 Min.   : 90.0   Min.   : 48.00   Min.   :26.00   Min.   :30.00  
 1st Qu.:121.2   1st Qu.: 75.00   1st Qu.:33.00   1st Qu.:39.00  
 Median :136.0   Median : 82.00   Median :37.00   Median :42.00  
 Mean   :137.2   Mean   : 83.36   Mean   :37.93   Mean   :43.05  
 3rd Qu.:148.0   3rd Qu.: 92.00   3rd Qu.:41.75   3rd Qu.:46.00  
 Max.   :250.0   Max.   :124.00   Max.   :56.00   Max.   :64.00  
    time.ppn             location      gender       frame    
 Min.   :   5.00   Buckingham:175   female:214   small : 98  
 1st Qu.:  93.75   Louisa    :191   male  :152   medium:172  
 Median : 240.00                                 large : 96  
 Mean   : 339.04                                             
 3rd Qu.: 480.00                                             
 Max.   :1560.00                                             
```


## Reordering rows

We now have a clean dataset to work with, congratulations!
Let's see how we can order the rows according to certain columns. We will use the function `arrange()` from the `dplyr` package.

### Sorting according to numerical column


```r
### order the rows by increasing age
dat.nona = dat.nona %>% arrange(age)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  193       77  49   3.9  4.31  19     61    119   118    70    32  38
2  146       79  41   3.6  4.76  19     60    135   108    58    33  40
3  230      112  64   3.6  4.53  20     67    159   100    90    31  39
4  193      106  63   3.1  6.35  20     68    274   165   110    49  58
5  170       69  64   2.7  4.39  20     64    161   108    70    37  40
6  164       71  63   2.6  4.51  20     72    145   108    78    29  36
  time.ppn   location gender  frame
1      300     Louisa female  small
2      240 Buckingham female medium
3     1440     Louisa   male medium
4       60 Buckingham female  small
5      120 Buckingham female medium
6     1080 Buckingham   male  small
```

If we want to order by decreasing age:


```r
dat.nona = dat.nona %>% arrange(desc(age))
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  165       94  69   2.4  4.98  92     62    217   160    82    51  51
2  301       90 118   2.6  4.28  89     61    115   218    90    31  41
3  226      279  52   4.3 10.07  84     60    192   144    88    41  48
4  227      105  44   5.2  5.71  83     59    125   150    90    35  40
5  240       88  49   4.9  4.92  82     63    170   180    86    41  46
6  271      121  40   6.8  4.57  81     64    158   146    76    36  43
  time.ppn   location gender  frame
1      180 Buckingham female  large
2      210     Louisa female medium
3      210     Louisa female  small
4      300     Louisa female medium
5      720 Buckingham female medium
6       10     Louisa female medium
```

### Sorting according to a column containing strings or factors


```r
dat.nona = dat.nona %>% arrange(gender)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  165       94  69   2.4  4.98  92     62    217   160    82    51  51
2  301       90 118   2.6  4.28  89     61    115   218    90    31  41
3  226      279  52   4.3 10.07  84     60    192   144    88    41  48
4  227      105  44   5.2  5.71  83     59    125   150    90    35  40
5  240       88  49   4.9  4.92  82     63    170   180    86    41  46
6  271      121  40   6.8  4.57  81     64    158   146    76    36  43
  time.ppn   location gender  frame
1      180 Buckingham female  large
2      210     Louisa female medium
3      210     Louisa female  small
4      300     Louisa female medium
5      720 Buckingham female medium
6       10     Louisa female medium
```

Given that there are many patients with the same value in this column, how can we order the patients within a certain category (for example, sorting the female patients by increasing age):


```r
dat.nona = dat.nona %>% arrange(gender, age)
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  193       77  49   3.9  4.31  19     61    119   118    70    32  38
2  146       79  41   3.6  4.76  19     60    135   108    58    33  40
3  193      106  63   3.1  6.35  20     68    274   165   110    49  58
4  170       69  64   2.7  4.39  20     64    161   108    70    37  40
5  149       77  49   3.0  4.50  20     62    115   105    82    31  37
6  226       97  70   3.2  3.88  20     64    114   122    64    31  39
  time.ppn   location gender  frame
1      300     Louisa female  small
2      240 Buckingham female medium
3       60 Buckingham female  small
4      120 Buckingham female medium
5      720 Buckingham female  small
6       90     Louisa female  small
```

> Order the rows by location, then gender, and decreasing weight!

<details>
<summary><b>Click for solution!</b></summary>


```r
dat.nona = dat.nona %>% arrange(location, gender, desc(weight))
head(dat.nona)
```

```
  chol stab.glu hdl ratio glyhb age height weight bp.1s bp.1d waist hip
1  192      109  44   4.4  4.86  43     64    325   141    79    53  62
2  235      109  59   4.0  7.48  62     63    290   175    80    55  62
3  203      299  43   4.7 12.74  38     69    288   136    83    48  55
4  193      106  63   3.1  6.35  20     68    274   165   110    49  58
5  180       84  69   2.6  5.20  40     68    264   142    98    43  54
6  228       92  37   6.2  4.64  58     61    256   190    92    49  57
  time.ppn   location gender  frame
1       60 Buckingham female  large
2      300 Buckingham female  large
3      240 Buckingham female  large
4       60 Buckingham female  small
5      240 Buckingham female medium
6      180 Buckingham female  large
```

</details>
