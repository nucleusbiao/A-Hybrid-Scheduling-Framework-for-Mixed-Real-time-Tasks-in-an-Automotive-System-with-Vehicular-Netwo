# A-Hybrid-Scheduling-Framework-for-Mixed-Real-time-Tasks-in-an-Automotive-System-with-Vehicular-Network


Preparation
------------------------------------------------
--Please visit https://www.mpa.ethz.ch/ to download real-time calculus toobox
--Read rtc/install.txt

	cd $MATLAB$/toolbox/rtc
	
	run rtc/rtc_install.m


Run code
------------------------------------------------------------------------------------------------------------
1. Test (Simple to use)

you can use USDPLUS_test.m, Exact_test.m, GA/testbyGA.m to get a scheduling result for a test quickly.
These files correspond to the algorithms in the article：UT-SD+， Exact， GA.
Because the codes of UT-SD and UT-SD+ are very similar, we only present the test codes of UT-SD+.
In these files, we define a set of random small task data for quick execution.

And the Online_test.m can be used to do a simple online assignment test.

1. The experiment code  (Load data or Create your data)

In data folder, there are the task data files that we used in the article. You can use in USD.m(UT-SD), 
USDPLUS.m(UT-SD+), Exact.m(Exact), and GA/GA.m(GA), and jobAssignOnline.m.

And you can use CreateTask.m or generate_task.mlx to generate your own data. Use them in the test code or the experiment code
