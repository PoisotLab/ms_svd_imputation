---
bibliography: [references.bib]
---

The «science of survival» concept has been introduced by Potter in 1970, and
consists in a reflection on the interdependence between a host and its
environment. Potter emphasizes the fact that the humanity's survival relies on
the surrounding biological and ecological systems, their proper functioning, and
can be influenced by their dynamics (Potter, 1970). The current pandemic of
SARS-CoV-2 is a stark reminder that movement of viruses through novel animal
hosts, and ultimately to human (zoonotic spillovers; Plowright *et al.*,
2017) requires that we understand the complexity of our biological surroundings.
Indeed, the fact that the majority of emerging infectious diseases are caused by
zoonotic pathogens from wildlife sources (Jones *et al.*, 2008) gives some
urgency to the task of predicting which viruses can be found in which hosts, so
as to provide guidance on where and what species to sample and where spillovers
are likely to happen(Johnson *et al.*, 2020, Albery *et al.*, 2020).

As seen with SARS-CoV and MERS-CoV epidemics, novel human infections by viruses
are representing a serious threat to global public health, and being able to
prevent future viral emergence now appears as a fundamental tool among our
society. Zoonotic dynamics usually involve three main stages: transmission
within the animal reservoir, cross-species spillover and transmission to human,
and finally, transmission among humans (Lloyd-Smith *et al.*, 2009). In the
past decades, substantial research effort has been put in studying and
predicting dynamics at the animal-human interface, but tracing back the ultimate
origin of novel zoonotic viruses remains a major difficulty (Becker *et*
al.}, 2020). Also, the main strategy adopted so far against infectious diseases
consists in taking actions after the emergence by increasing the health
infrastructures and vigilance, as well as developing vaccines or medical
treatments (Han *et al.*, 2016).  

As suggested by Han *et el.*, a more efficient approach would be
anticipatory. Yet an anticipatory approach can be limited by lack of suitable
data, and as Becker *et al.* (2020) highlighted, by disagreement between
models. The task of predicting possible host-virus interactions would therefore
benefit from adding methods that allow imputation, and can produce results that
are easily added to ensemble models. Here, we explore an approach focusing on
the first stage of zoonoses dynamics, by using the Singular Value Decomposition
(SVD) as an imputation method for identifying unobserved host-virus
interactions, acting as potential intermediates hosts in diseases transmissions.

# Theory

Singular Value Decomposition (SVD) is a linear algebra technique used to
decompose a data matrix in a product of three matrices (Ginebreda *et al.*,
2019):

$$\mathbf{X} =  \mathbf{U \Sigma V}^T$$ {#eq:svd}

Where $\mathbf{X}$ is a $mxn$ data matrix ($m \ge n$), $\mathbf{U}$ is
an unitary $mxm$ matrix containing the left singular vectors, $\mathbf{V}$
is an unitary $nxn$ matrix containing the right singular vectors and
$\mathbf{\Sigma} $ is a diagonal matrix containing the singular values ordered
in decreasing order of importance, in regard of the quantity of information that
they present.

This process allows data reduction by finding key correlations among entries and
then by approximating the original matrix.  

\paragraph{Optimal Truncation} of the SVD at rank $r$  of the singular values
will allow data reduction while keeping enough information to obtain a balance
between complexity and accuracy within the model. We ran the analyses in
*Julia* 1.5.1 (Bezanson *et al.* 2017). Truncation at rank $r$ was
performed by setting values $\mathbf{\Sigma}_{(r+1)..m}$ to 0 (we note the
resulting vector $^{(r)}\mathbf{\Sigma}$), and the resulting low-rank
approximation was obtained by

$$^{(r)}\mathbf{X} =  \mathbf{U} \, ^{(r)}\mathbf{\Sigma \, V}^T$$ {#eq:lowrank}

# Method

## Model structure

For each non-interactions in the dataset (see next section), the model then
assigns an initial value to it and performs iteratively the SVD at chosen rank,
until it reaches convergence. During this step, the cells in the matrix that are
*not* being imputed are kept at the actual value. We capped the maximal
number of iterations at 50, even though the value of the imputed cells stopped
changing (defined as a step-wise change lower than $10\times \epsilon$) after
less than 10 steps in most cases. The initial value we picked for this
illustration is the connectance of the global host-virus interaction dataset,
which amounts to the probability that any pair of organisms are found to
interact (0.03). Other schemes to impute the initial value are possible, for
example by relying on the relative degree of both species, but this would
require more guesswork or more assumptions, in addition to possibly being
sensitive to biases in the preferential sampling of some groups.

## Dataset

The method has been applied to the dataset assembled by the VERENA consortium on
mammalian viruses
(`https://github.com/viralemergence/virionette/blob/master/03_interaction_data/virionette.csv).
This dataset regroups interactions between 710 mammalian hosts and 72 viruses,
and the aggregation process has been described in Becker *et al.* (2020).
Specific attention has been paid to betacoronaviruses, a viral genus at high
risk of spillover, and to their potential bat hosts, a mammalian order known to
be evolutionary implied in the main viruses zoonotic historical epidemics (Ren
*et al.*, 2006).

## Method validation

For every rank of the SVD, we examined the top 10 non-interactions with the
highest score (even though the model technically returns a score for all
non-interactions). Since there have recently been a number of empirical
investigations regarding host-virus associations, novel bat hosts of
betacoronaviruses have been described following the assembly of the Becker
*et al.* (2020) dataset used for the model. Therefore, we were able to
examine and compare which of these novel identified bat hosts for
betacoronaviruses have also been recovered by SVD, at all ranks between 1 and 5,
to confirm the accuracy of this technique.

# Results and Discussion

First, we report the top 10 likely hosts for betacoronaviruses, which are ranked
by their final value post imputation; larger values should indicate that the
interactions are more likely to be possible. We report the novel hosts
(identified post Becker *et al.* 2020, according to
`https://www.viralemergence.org/betacov`). These results are presented in Table
1 - the novel hosts are presented in **bold**. Using a rank 2 approximation of
the dataset, we have 5 novel hosts, and 4 identified as "suspected" hosts by the
Becker *et al.* (2020) ensemble model, currently lacking empirical evidence.
This suggests that rank 2 contains the most information about the processes
generating the data, and can therefore be used to infer other associations.

Based on this information, we have extracted the 10 highest scoring interactions
across the entire matrix at rank 2 (Table 2). The results demonstrates that
within the entire dataset, including all mammalian hosts and viruses' genus, 5
out of the 10 highest scoring interactions are involving bat hosts (presented in
*italic*), and 8 out of the 10 interactions are involving the lyssavirus
genus. This genus includes the rabies virus (RABV), and other neurotropic
rabies-related viruses (Warrell *et al.*, 2004).

Once those results were obtain, further investigations in the form of literature
surveys allowed to identify that the interaction between *Pipistrellus*
abramus} and lyssaviruses has already been noted by Hu *et al.* in 2018;
Shipley *et al.* (2019) reported lyssavirus prevalence in the genus
*Pipstrellus*, *Myotis*, and *Rhinolophus*. Other confirmed
hosts of lyssaviruses are *Sus scrofa* (Satou *et al.*, 2004), and
*Rattus norvegicus* (*e.g.* Wang *et al.*, 2014).
Surveillance for novel lyssaviruses infections is of great public health
interest, since the rabies virus is fatal in all cases, once the onset of
clinical symptoms has started (Banyard *et al.*, 2017). Although it is
recognized that bats are identified as reservoir hosts for lyssaviruses, the
mechanism allowing the maintenance of the virus in those populations is still
poorly understood (Banyard *et al.*, 2017), and these predictions of
interactions might serve as guidance in the monitoring of new infections. 

The two non-lyssaviruses associations have been previously reported in the
literature (*Sus scrofa* and orbivirus by Belaganhalli *et al.*,
2015; *Capra hircus* and the equine encephalomyelitis caused by an
alphavirus as early as Purcell *et al.*, 1976). This suggests that
Singular Value Decomposition of available data on host-virus associations can
uncover results that have been reported in the primary literature, but not
incorporated in the main databases used in the field; based on the fact that the
majority of the top 10 overall associations were able to be validated from the
literature, we suggest that interactions that have no empirical evidence could
be targets for additional sampling.

| Rank 1                   | Rank 2                    |
|--------------------------|---------------------------|
| **Artibeus jamaicensis** | **Hipposideros pomona**   |
| **Scotophilus kuhlii**   | **Scotophilus kuhlii**    |
| Molossus rufus           | **Artibeus jamaicensis**  |
| Sturnira lilium          | Carollia brevicauda       |
| **Desmodus rotundus**    | Chaerephon pumilus        |
| Glossophaga soricina     | Molossus rufus            |
| Eptesicus fuscus         | Glossophaga soricina      |
| Tadarida brasiliensis    | **Desmodus rotundus**     |
| Myotis nigricans         | Sturnira lilium           |
| Myotis lucifugus         | **Hipposideros larvatus** |

| Hosts species          | Viruses genus |
|------------------------|---------------|
| Sus scrofa             | Lyssavirus    |
| *Hipposideros armiger* | Lyssavirus    |
| Rattus norvegicus      | Lyssavirus    |
| Myodes glareolus       | Lyssavirus    |
| *Pipistrellus abramus* | Lyssavirus    |
| Sus scrofa             | Orbivirus     |
| Capra hircus           | Alphavirus    |
| *Rhinolophus sinicus*  | Lyssavirus    |
| *Myotis ricketti*      | Lyssavirus    |
| *Rhinolophus affinis*  | Lyssavirus    |

# Conclusion and future work

Being able to identify intermediate animal hosts for potential zoonotic
pathogens is an important step in the fight against potential threats to global
public health. Using SVD as an imputation method to predict those interactions
has demonstrated its potential to achieve this goal by correctly identifying the
majority of the most likely associations, as validated by literature surveys,
and by suggesting interactions with no empirical evidence as targets for
additional sampling. Host-virus associations are a challenging imputation
problem, because organized datasets are scarce -- as a result, a lot of missing
associations are reported in the literature, but not available in an easily
usable format. Yet this also presents an opportunity to validate the performance
of recommender systems that is far more interesting than cross-fold or
leave-one-out validation: the existence of these interactions in the literature
can provide validation on data that have never been used in the modeling
process, and therefore provide an accurate estimate of how frequently existing
interactions are identified. By this measure, that most of the top 10
recommendations on this dataset were validated through *de novo* sampling
(for bat hosts of betacoronaviruses) or by a literature survey (for the global
dataset) is a strong indication that SVD is able to uncover likely host-virus
pairs.

Future work on the use of SVD for virus host associations will have to adress
the question of the initial value used in the imputation process. As of now, we
relied on the average number of interactions in the matrix; yet this can
overestimate the importance of viruses with a narrow host range, or
underestimate the importance of generalist viruses. For this reason, we are
confident that the performance of the approach can further be improved by
fine-tuning the initial value used for imputation. A promising avenue in this
regard is to rely on Stock *et al.* (2017) work on linear filtering --
this provides a convenient way to assign weights to various aspects of network
structure, and has been revealed to provide a good baseline estimate of how
likely it is that a missing interaction actually exists. Combining an accurate
model for the initial value with the SVD imputation is likely to generate
predicted interactions that are strong candidates for empirical validation.

# References

Albery, G. F., Eskew, E. A., Ross, N., \& Olival, K. J. (2020). Predicting the global mammalian viral sharing network using phylogeography. Nature communications, 11(1), 1-9.

Banyard, A. C., \& Fooks, A. R. (2017). The impact of novel lyssavirus discovery. Microbiology Australia, 38(1), 17-21.

Becker, D., Albery, G. F., Sjodin, A. R., Poisot, T., Dallas, T., Eskew, E. A., \& Carlson, C. J. (2020). Predicting wildlife hosts of betacoronaviruses for SARS-CoV-2 sampling prioritization. bioRxiv.

Belaganahalli, M. N., Maan, S., Maan, N. S., Brownlie, J., Tesh, R., Attoui, H., \& Mertens, P. P. (2015). Genetic characterization of the tick-borne orbiviruses. Viruses, 7(5), 2185-2209.

Bezanson, J., Edelman, A., Karpinski, S., \& Shah, V. B. (2017). Julia: A fresh approach to numerical computing. SIAM review, 59(1), 65-98.

Ginebreda, A., Sabater-Liesa, L., \& Barceló, D. (2019). Quantification of ecological complexity and resilience from multivariate biological metrics datasets using singular value decomposition entropy. MethodsX, 6, 1668-1676.

Han, B. A., \& Drake, J. M. (2016). Future directions in analytics for infectious disease intelligence: toward an integrated warning system for emerging pathogens. EMBO reports, 17(6), 785-789.

Hu, S. C., Hsu, C. L., Lee, M. S., Tu, Y. C., Chang, J. C., Wu, C. H., ... \& Tu, W. J. (2018). Lyssavirus in japanese pipistrelle, taiwan. Emerging infectious diseases, 24(4), 782.

Johnson, C. K., Hitchens, P. L., Pandit, P. S., Rushmore, J., Evans, T. S., Young, C. C., \& Doyle, M. M. (2020). Global shifts in mammalian population trends reveal key predictors of virus spillover risk. Proceedings of the Royal Society B, 287(1924), 20192736.

Jones, K. E., Patel, N. G., Levy, M. A., Storeygard, A., Balk, D., Gittleman, J. L., \& Daszak, P. (2008). Global trends in emerging infectious diseases. Nature, 451(7181), 990-993.

Lloyd-Smith, J. O., George, D., Pepin, K. M., Pitzer, V. E., Pulliam, J. R., Dobson, A. P., \& Grenfell, B. T. (2009). Epidemic dynamics at the human-animal interface. science, 326(5958), 1362-1367.

Plowright, R. K., Parrish, C. R., McCallum, H., Hudson, P. J., Ko, A. I., Graham, A. L., \& Lloyd-Smith, J. O. (2017). Pathways to zoonotic spillover. Nature Reviews Microbiology, 15(8), 502-510.

Potter, V. R. (1970). Bioethics, the science of survival. Perspectives in biology and medicine, 14(1), 127-153.

Pursell AR, Peckham JC, Cole JR, Stewart WC, Mitchell FE, (1972). Naturally occurring and artificially induced eastern encephalomyelitis in pigs. Journal of the American Veterinary Medical Association, 161(10):1143-1147.

Ren, W., Li, W., Yu, M., Hao, P., Zhang, Y., Zhou, P., ... \& Wang, L. F. (2006). Full-length genome sequences of two SARS-like coronaviruses in horseshoe bats and genetic variation analysis. Journal of General Virology, 87(11), 3355-3359.

Sato, G., Itou, T., Shoji, Y., Miura, Y., Mikami, T., Ito, M., ... \& Ito, F. H. (2004). Genetic and phylogenetic analysis of glycoprotein of rabies virus isolated from several species in Brazil. Journal of Veterinary Medical Science, 66(7), 747-753.

Shipley, R., Wright, E., Selden, D., Wu, G., Aegerter, J., Fooks, A. R., \& Banyard, A. C. (2019). Bats and viruses: Emergence of novel lyssaviruses and association of bats with viral zoonoses in the eu. Tropical medicine and infectious disease, 4(1), 31.

Stock, M., Poisot, T., Waegeman, W., \& De Baets, B. (2017). Linear filtering reveals false negatives in species interaction data. Scientific reports, 7(1), 1-8.

Wang, L., Tang, Q., \& Liang, G. (2014). Rabies and rabies virus in wildlife in mainland China, 1990–2013. International Journal of Infectious Diseases, 25, 122-129.

Warrell, M. J., \& Warrell, D. A. (2004). Rabies and other lyssavirus diseases. The Lancet, 363(9413), 959-969.
