#Bayesian Compressive Sensing#
This set of Matlab (7.0) functions contain the core code for generating 
Fig.2 and Figs.4(a,b) of the following paper:

"Bayesian Compressive Sensing" 
Shihao Ji, Ya Xue, and Lawrence Carin, IEEE Trans. Signal Processing, vol. 56, no. 6, June 2008.

and Fig.2 and Fig.3 of the following paper:

"Multi-Task Compressive Sensing" 
Shihao Ji, David Dunson, and Lawrence Carin, IEEE Trans. Signal Processing (Accepted), 2008.


##BCS_demo##

1) Fig2.m ---> generate Fig.2
     The following two Matlab files from l1-magic are required for BP implementation:
     l1qc_logbarrier.m
     l1qc_newton.m
     which can be found at: http://www.acm.caltech.edu/l1magic/

2) Fig4_ab.m ---> generate Figs.4(a,b)
     multi_random_measures.m     ----> generate the "Random" curve for Fig.4a
     multi_optimized_measures.m  ----> generate the "Optimized" curve for Fig.4a
     multi_approx_measures.m     ----> generate the "Approx." curve for Fig.4a


##MT_CS_demo##

***A bug was fixed on Aug. 03, 2008 in MT_CS.m for the cases where signals are dramatic undersampled.***
1) Fig2.m ---> generate Fig.2

2) Fig3.m ---> generate Fig.3
     multi_runs_75.m     ----> generate the "MT 75%" curve for Fig.3
     multi_runs_50.m     ----> generate the "MT 50%" curve for Fig.3
     multi_runs_25.m     ----> generate the "MT 25%" curve for Fig.3


The code is still in development stage. If you have any comments or bug reports, 
you are welcome to contact Shihao Ji at shji@ece.duke.edu.

Please check back at http://www.ece.duke.edu/~shji/BCS.html for updates.

Shihao Ji
Duke University, July 8, 2007

#LICENSE#
Distribution and use of this code is subject to the following agreement:

This Program is provided by Duke University and the authors as a service
to the research community. It is provided without cost or restrictions, 
except for the User's acknowledgement that the Program is provided on an 
"As Is" basis and User understands that Duke University and the authors 
make no express or implied warranty of any kind.  Duke University and the
authors specifically disclaim any implied warranty or merchantability or 
fitness for a particular purpose, and make no representations or warranties 
that the Program will not infringe the intellectual property rights of 
others. The User agrees to indemnify and hold harmless Duke University and
the authors from and against any and all liability arising out of User's 
use of the Program.
