***---------------------------------------------***
***												***
*** Prevalence and Trend Estimation from 		***
*** Observational Data with Highly Variable 	***
*** Post-Stratification Weights					***
***												***
*** Y. Vandendijck, C. Faes, N. Hens			***
***												***
*** Supplementary Materials: Appendix C			***
***												***
***---------------------------------------------***

*** Fit the model;
***--------------;

* Read data;
proc import datafile='F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\PrevalenceExample.txt'
     out=dataInput
     dbms=dlm
	 replace;
     datarow=2;
	 *GUESSINGROWS=18;
run;
data dataInput;
	set dataInput;
	hcopy = xh;
run;

* LIN - model;
proc glimmix data=dataInput method=RSPL;
	title 'LIN model';
	parms (0.02);
	class hcopy;
	model resp/nh = xh /solution dist=binomial;
	random intercept /subject=hcopy solution;
	output out=DataOutput1 predicted=eta predicted(ILINK)=prob variance=var;
	ods output CovParms=DataOutput2;
run;
proc export data= DataOutput1 
            outfile= 'F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\InputSAS1.txt' 
            dbms=TAB replace;
     putnames=yes;
RUN;
proc export data= DataOutput2 
            outfile= 'F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\InputSAS2.txt' 
            dbms=TAB replace;
     putnames=yes;
RUN;


* XRE - model;
proc glimmix data=dataInput method=RSPL;
	title 'XRE model';
	parms (0.02);
	class hcopy;
	model resp/nh= /solution dist=binomial ;
	random intercept /subject=hcopy solution;
	output out=DataOutput1 predicted=eta predicted(ILINK)=prob variance=var;
	ods output CovParms=DataOutput2;
run;


* NPAR - model;
proc glimmix data=dataInput method=RSPL;
	title 'NPAR model';
	parms (0.001) (0.02);
	class hcopy;
	model resp/nh = xh /solution dist=binomial;
	random xh /type=rsmooth solution knotmethod=equal(18) knotinfo;
	random intercept /subject=hcopy solution;
	output out=DataOutput1 predicted=mu predicted(ILINK)=prob variance=var;
	ods output CovParms=DataOutput2;
run;






*** Bootstrap method;
***-----------------;

* Read data;
proc import datafile='F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\BootstrapData.txt'
     out=dataInput
     dbms=dlm
	 replace;
     datarow=2;
	 GUESSINGROWS=18;
run;
data dataInput;
	set dataInput;
	hcopy = xh;
	group=1;
run;

* LIN - model;
proc glimmix data=dataInput method=RSPL;
	by b;
	title 'LIN model';
	parms (0.02);
	class hcopy;
	model y/nh = xh /solution dist=binomial;
	random intercept /subject=hcopy solution;
	output out=DataOutput1 predicted(ILINK)=prob;
	nloptions maxtime=10;
run;
data DataOutput1;
	set DataOutput1;
	if prob=. then prob=999;
run;
proc export data= DataOutput1 
            outfile= 'F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\BootstrapInputSAS1.txt' 
            dbms=TAB replace;
     putnames=yes;
RUN;


*** Jackknife method;
***-----------------;

* Read data;
proc import datafile='F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\JackknifeData.txt'
     out=dataInput
     dbms=dlm
	 replace;
     datarow=2;
	 GUESSINGROWS=18;
run;
data dataInput;
	set dataInput;
	hcopy = xh;
	group=1;
run;

* LIN - model;
proc glimmix data=dataInput method=RSPL;
	by b;
	title 'LIN model';
	parms (0.02);
	class hcopy;
	model y/nh = xh /solution dist=binomial;
	random intercept /subject=hcopy solution;
	output out=DataOutput1 predicted(ILINK)=prob;
	nloptions maxtime=10;
run;
data DataOutput1;
	set DataOutput1;
	if prob=. then prob=999;
run;
proc export data= DataOutput1 
            outfile= 'F:\PhD\Weighting\Paper\Paper major revision 1 files\Appendix C\JackknifeInputSAS1.txt' 
            dbms=TAB replace;
     putnames=yes;
RUN;
