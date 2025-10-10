# Epistatic Network and Dominance Prediction in Seasonal Influenza A (H3N2)

## Overview

Our study aimed to evaluate how the **epistatic network among eight gene segments** shapes the **evolutionary dynamics of the haemagglutinin (HA) gene** in seasonal influenza A (H3N2).  
We quantified **synonymous and nonsynonymous genetic distances** from vaccine strains and developed predictive models to estimate the **dominance of emerging strains** across vaccine periods.

---

## [Bayesian Renaissance Counting (RNSC)]

![Figure 1: Bayesian Renaissance Counting (RNSC)](3_images/githubFig1.jpg)

We applied **Bayesian Renaissance counting (RNSC)** within **BEAST** to estimate synonymous and nonsynonymous substitution distances in a time-resolved phylogeny.  
These distances quantitatively represent the **degree of genetic heterogeneity** between circulating variants and the corresponding vaccine strain, reflecting **immune-driven divergence**.

In the dominance prediction model, strains are visually outlined by their dominance status (see figure).  
**Dominance** was defined by the strain’s **existence before and persistence after a specific vaccine period**, indicating successful adaptation and continued transmission.

---

## [Training/Testing Datasets for Mixed-Effect Logistic Regression]

![Figure 2: Mixed-Effect Logistic Regression Framework](3_images/githubFig2.jpg)

We constructed a **mixed-effect logistic regression model** to compare the predictive power of genetic information derived from the **HA gene alone** versus the **multi-segment epistatic network**.  
This regression-based framework quantitatively estimated the **likelihood of dominance** for H3N2 strains across consecutive vaccine updates.

Importantly, **adaptive mutations in the NA and NS segments**, which are likely associated with **immune escape**, significantly enhanced the **predictive accuracy for HA dominance**.

---

## [Training/Validation/Testing Datasets for Machine Learning Models]

![Figure 3: Machine Learning and Deep Learning Comparison](3_images/githubFig3.jpg)

To extend and validate the framework, we compared the logistic mixed-effects baseline model with **machine learning** approaches (Random Forest, SVM, LASSO, Ridge Regression) and **deep learning** (DNN).  
These models further improved the robustness of dominance prediction by incorporating **multi-segment genetic signals** and capturing **nonlinear epistatic interactions**.

---

## Repository Structure

