hold on
load('solution_values5.mat')
load('solution_values4.mat')
plot(bpp_solution, psnr_solution, '--*'); % reference curve Chapter 5
plot(bpp_solution_ch4, psnr_solution_ch4, '--+'); % reference curve Chapter 4