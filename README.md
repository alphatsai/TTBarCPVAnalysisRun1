# [Regression] Top-qaurk CP violation in LHC data

<div style="text-align: center;" markdown="1"><a href="http://dx.doi.org/10.1007/JHEP03(2017)101">JHEP 03 (2017) 101</a></div>

*"Search for CP violation in tt production and decay in proton-proton collisions at √s= 8 TeV"*

*<div style="text-align: center;" markdown="1"> `minimum chi-square` `maximum likelihood` `big data` `statistics`</div>*

Original website: https://hackmd.io/s/BJl6iAVVz

## Introduction
One of popular question about our universe is the large asymmetric amounts between ***matter*** and ***antimatter***, while it violates the common knowledge that the universe is balanced and symmetric in any physics mechanism, e.g. the traveling photon can tranform and split (decay) to electron ($e^-$) and anti-electron ($e^+$), called *positron* as well, and the total charge, spin and momentum are conserved. However, the observation beyonds the expectation, the most the anti-matter is gone, and matter (you, me, your dog and everything around us) is well staying in this universe, i.e. we are not omitted by the antimatter which sounds like the movie [*"Da Vinci Code - Angels & Demons (2009)"*](http://www.imdb.com/title/tt0808151/) (many physicists are killed at CERN :cry:)

The asymmetry is called the ***charge-parity (CP) violation*** , where the parity is about the mirror symmetry of spin and momentum. The ***Standard Model (SM)*** is the possible physics model to explain the phenomena, it predicts the nature CP volation can be observed from the decay of fundamental particles. Many experiments prove the model, but the amount is still far away from the observation in the universe. 

The top-quark's CP violation has not been observed, since the SM predicts it is too small to be asymmetric. Thus, this is the first search with top-quark via the proton-proton collision produces the top-quark pair, $\text{pp}\to t\bar{t}$. If any asymmetry is observed, this may be a light of the explaination.

The measurment is a simple data counting in terms of the observables $O_i$:

$$
\begin{equation}
    A_\mathrm{CP}(O_{i})=\frac{N_{O_{i}>0}-N_{O_{i}<0}}{N_{O_{i}>0}+N_{O_{i}<0}}\ ,
\end{equation}
$$

which is called ***asymmetric CP***; $N_{O_{i}>0}$ is number of data has positive value of $O_i$, while $N_{O_{i}<0}$ is for the negative value. As the $A_{\text{CP}}\neq0$, the measurement disagrees with the SM zero prediction, and the large asymmetry is possible from other mechanism. According to an alternative theory, there are four independent observables consisted of kinematics of the finals decays of the top-quark pair prodction in LHC big data, i.e. 2 b-quark jets ($b$, $\bar{b}$), 2 non-b quark jets ($j_1,j_2$), one electron/muon ($e, \mu\in\ell$) and neutrino ($\nu$).

<div style="text-align: center;" markdown="1"><img src="https://i.imgur.com/cRifnfI.png" height="250" ></div>
</br>

The observable are defined as:

$$
\begin{split}  
&O_{2}= (\vec{p}_{b}+\vec{p}_{\bar{b}})\cdot(\vec{p}_{\ell}\times\vec{p}_{j_{1}}),\\
&O_{3}= Q_{\ell}\,\vec{p}_{b}\cdot(\vec{p}_{\ell}\times\vec{p}_{j_{1}}),\\
&O_{4}= Q_{\ell}\,(\vec{p}_{b}-\vec{p}_{\bar{b}})\cdot(\vec{p}_{\ell}\times\vec{p}_{j_{1}}),\\
&O_{7}= (\vec{p}_{b}-\vec{p}_{\bar{b}})_{z}(\vec{p}_{b}\times\vec{p}_{\bar{b}})_{z}\ ,
\end{split}
$$

where $\vec{p}$ are the three-momenta of the final-state particles; the subscript $z$ indicates a projection along the direction of the counterclockwise rotating proton beam, defined to be the $+z$ direction in the CMS coordinate system; and $Q_{\ell}$ is the electric charge of $\ell$. Note that the sign of the observable is the only information needed to measure $A_\text{CP}$.

The irreducible resolution of detector and particle reconstrustion algorithm dilutes the measuremnts. Moreover, the sensitivity of the observables depends on whether distinguishable objects are involved in their definition. For instance, the $b$ quark jet charges need to be distinguished for $O_3$ and $O_4$, but not for $O_2$ and $O_7$.



## Techniques
This is the new measurement with LHC data, the techniques for extracting the interesting parameters, i.e. $A_{\text{CP}}$, are well designed. Since the measurement is only for the particular physics process, the techniques to deal with the background noise includes the [**Data selection**](#1-Data-selection),  [**Minimum $\chi^2$ method**](#2-Minimum-chi2-method) and [**Maximum likelihood regression**](#3-Maximum-likelihood-regression). The first filters the most background and defines the ***"data-driven control sample"*** for background's fake $A_{\text{CP}}$. The second the *statistical method* for extracting the feature hiding in the data with the physics assumption. At the third, the regression method is used to estimate the signal and backgroud noise from selected data with a particular *likelihood function*. In the end, with subtracting the predicted background, the $A_{\text{CP}}$ can be extracted from data. 

### 1. Data selection
As the figured illustrated in [Introduction](#Introduction) and the listed observables, the used objects (reconstructed particles) are 2 b-quark jets, the highest non-b quark jet and an isolated electron/muon, which luckily kills many background noise in LHC data. The selected data still contains a few irreducible background noise. However, since the background is expected to have zero $A_{\text{CP}}$ due to the random combination, we extract the background kinematic distributions (shapes) by using the outlier data which has no any b-quark jets, called ***data-driven*** background. It is validated with ***Monte Carlo (MC)*** samples. The background subtraction from selected data is done with [Maximum likelihood regression](#3-Maximum-likelihood-regression) to obtained the expected number of background yields. Here is one of kinematic distribution which is going to be used for signal and background prediction:

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/xhBKmvj.png" height="300" >](https://i.imgur.com/xhBKmvj.png)</div>

The "events" means the data quantity which is used in HEP to represent a physics interation happened; The colorful filled histograms are made by MC samples.

### 2. Minimum $\chi^2$ method
The observable is sensitive to whether objects are distiguishable without random combination. The challenge in our case is to distinguish the $b$- and $\bar{b}$-quark jets, since there is only the classification algorithm for b-qaurk and other quarks instead of its charge.  

$$
 \chi^2 = \left(\frac{m_{bjj}-m_{t}}{\sigma_{t}}\right)^2 + \left(\frac{m_{jj}-m_{W}}{\sigma_{W}}\right)^2 \ ,
$$

where $m_{bjj}$ ($m_{jj}$) is the invariant mass of a b-tagged jet and two non-b-tagged jets (two non-b-tagged jets); $m_{t}$ ($m_{W}$) and $\sigma_{t}$ ($\sigma_{W}$) is the mean value and width of top quark ($W$ bosson) mass which are the well-known constants. By scaning the all objects, the combination of $m_{bjj}$ and $m_{jj}$ having the minimum $\chi^2$ is the proper choise. In the meanwhile, the charge of b-quark jets can be distiguished as the figure illustracted in [Introduction](#Introduction). The method is validated with MC samples, the efficiency and wrong-tagging rate are shown as 

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/JzkTmQr.png" height="300" >](https://i.imgur.com/JzkTmQr.png)</div>

We choose $\chi^2_{min}=40$ for selecting the data having good $t\bar{t}$ pair.  

### 3. Maximum likelihood regression
To estimate the signal and background yields in selected data, we use the likelihood function to do the regression fit. The probability distribution function of signal and background are obtained by MC and data-driven, respectively.  


#### 3.1. Defination
The case for estimating the data yield, the normal ***likelihood function***, which multiples the all the probabilities in terms of variable $x$ for all data, is multiplied the ***Possion distribution function***, since the sum of expected singal and background yields, $\mu=n_s+n_b$, is expected to be deviated from observed total yields, $N$, follow the *Poisson probability*. It is  defined as

$$
\begin{split}
L&=\frac{\mu^Ne^{-\mu}}{N!}\prod_{i=1}^N P(x_i|\mu) \\
&=\frac{(n_s+n_b)^Ne^{-(n_s+n_b)}}{N!}\prod_{i=1}^N{p(n_s)p(x_i|s)+p(n_b)p(x_i|b)} \\
&=\frac{e^{-(n_s+n_b)}}{N!}\prod_{i=1}^N{n_sp(x_i|s)+n_bp(x_i|b)}\ ,
\end{split}
$$

where $n_s$ and $n_b$ are the target parameters; $p(x|s,b)$ is the probability distribution of $x$ of the signal (s) or background (b) data in selected region. In the computing case, the optimization is usually looking for the minimum instead of the maximum. Thus we use logarithm of the likelihood to have the minimum solution: 

$$
S=-\ln{L} = n_s+n_b+\sum_{i=1}^N\ln\left(n_sp(x_i|s)+n_bp(x_i|b)\right)-\ln{N!}\ ,
$$

where the $-\ln{N!}$ can be ignored in the fitting, since it is a constant.

#### 3.2. Validation
The variable $x$ is used the distribution of *invariant mass of the $e$/$\mu$ and $b$*, denoted $m_{\ell b}$. The observed events (data) and the mean of predicted (fitted) events with $\pm1\sigma$ statistical (second term) and systematic (third term) uncertainties are as following list:

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/RVskFzA.png" height="100" >](https://i.imgur.com/RVskFzA.png)</div>

The validation is comparing the data and fitted distribution by ***Goodness of fit*** method. The *goodness of fit* method is caculating the $\chi^2/d.o.f.$ for each data with respect to the fitted distribution. The deviation of each data is assumed with following the ***Gaussian distribution***. 

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/zPxpleW.png" height="350" >](https://i.imgur.com/zPxpleW.png)</div>

## Results

### 1. Sensitivity
The measurements are used the reconstructed objects, thier resolutions are limited by detector and algorithms. The results is possible to be corrected to origainal $A_{\text{CP}}$ by ***"dilution factor"*** . We estimates the factors with *MC samples* by comparing the truth and reconstruction results. As following figure, the ***gradient*** of diluted asymmetric, called ***"effective $A_{\text{CP}}$"***  denoted $A'_{\text{CP}}$, is the scalar for correct the observation. 

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/GfJtnze.png" height="350" >](https://i.imgur.com/GfJtnze.png)</div>

### 2. Effective $A_\text{cp}$
The final $A'_{\text{CP}}$ and corrected $A_{\text{CP}}$ are covering the zero within the uncertainties, i.e. there is no new physics beyond the SM prediction... :cry:. At least this new measurement open another direction to search the unsolved puzzle. The measurement may have different behavior in much higher energy with improved techniques.

**<div style="text-align: center;" markdown="1">**[<img src="https://i.imgur.com/PkAIbLi.png" height="350" >](https://i.imgur.com/PkAIbLi.png)</div>



## References
- Physics dissertation, [dx.doi.org/10.6342/NTU201700205](http://www.airitilibrary.com/Publication/alDetailedMesh1?DocID=U0001-2301201711464200)
- [JHEP 03 (2017) 101](http://dx.doi.org/10.1007/JHEP03(2017)101)
- Github : https://github.com/juifa-tsai/TTBarCPVAnalysisRun1

