# Two-level optimization to estimate muscle and knee contact forces

Code containing a two-level optimization algorithm to estimate musculoskeletal parameters consistent with knee contact forces. This is part of the study published in the following article:

[Serrancolí, G; Kinney, A; Fregly, B.J. Influence of Musculoskeletal Model Parameter Values on Prediction of Accurate Knee Contact Forces during Walking. Medical Engineering & Physics 2020. 85: 35-47.](https://www.sciencedirect.com/science/article/pii/S1350453320301375?via%3Dihub "Reference article")

The main function to run the two-level optimization is [mainNestedOptimization_3trials.m](mainNestedOptimization_3trials.m), which allows you to estimate either optimal fiber lengths, tendon slack lengths, moment arms or a combination of all three. The user can reproduce the results of the mentioned article (9 two-level optimizations) by changing the options variables: Options.ma_opt, Options_lM0 and Optons.lTs.

If the user wants to run only the inner level optimization, one can run the function [mainInnerLevelOpt.m](mainInnerLevelOpt.m). The user needs to select the parent calibration (Options.parentCalibration), i.e. the number of the two-level optimization used to obtain the musculoskeletal parameters, and the walking cycle to estimate the muscle and knee contact forces (Options.patterncycle).




