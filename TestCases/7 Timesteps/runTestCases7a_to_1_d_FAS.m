diary 'TestCase_7a_to_te_FAS.txt'

%%
fprintf('\n FASTestCase_7a_100steps_24_1to5_V_ho \n');
for i=1:5
   run FASTestCase_7a_100steps_24_1to5_V_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_7b_80steps_24_1to5_V_ho \n');
for i=1:5
   run FASTestCase_7b_80steps_24_1to5_V_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
fprintf('\n FASTestCase_7c_60steps_24_1to5_V_ho \n');
for i=1:5
   run FASTestCase_7c_60steps_24_1to5_V_ho
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_7d_40steps_24_1to5_V_ho \n');
for i=1:5
   run FASTestCase_7d_40steps_24_1to5_V_ho
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

fprintf('\n FASTestCase_7e_20steps_24_1to5_V_ho \n');
for i=1:5
   run FASTestCase_7e_20steps_24_1to5_V_ho
end

%%
figure; % Lazy way of getting notice of the end of the simulations
diary off
