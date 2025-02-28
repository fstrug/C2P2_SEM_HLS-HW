#include <iostream>
#include <math.h>

// Unoptimized
#define SIZE 200
extern "C"{
void lin_fit(const int *x, const int *y, const int *y_err, const int *last, const int N){
    
    std::cout << "Entries: " << N << std::endl;
    // Loop through entries and calculate sums
    float S_x = 0;
    float S_y = 0;
    float S = 0;
    int x_i = 0;
    int y_i = 0;
    float y_err_i = 0;
    float y_err_sq = 0;
    int event_start = 0;
    int event_end = 0;
    float a = 0;
    float b = 0;
    
    
    float sigma_a2_partial, sigma_a2, sigma_b2;
    float chi = 0;
    float chi_2 = 0;
    float chi2_ndof = 0;
    
     // The user-provided values for the LOOP_TRIPCOUNT directive are used only for reporting,
    for (int i=0; i < N + 1; i++){
        #pragma HLS LOOP_TRIPCOUNT min=SIZE max=SIZE avg=SIZE
        if (i == N){
            break;
        }
        // Read in values to memory
        x_i = x[i];
        y_i = y[i];
        y_err_i = y_err[i];
        y_err_sq = y_err_i * y_err_i;
        // Calculate
        S_x += x_i/(y_err_sq);
        S_y += y_i/(y_err_sq);
        S += 1.0 / (y_err_sq);
        
        if (last[i] == 1){
            event_end = i;
            float S_tt = 0;
            float t_i = 0;
            float b_partial = 0;
            for (int j=event_start; j <= event_end; j++){
                #pragma HLS LOOP_TRIPCOUNT min=SIZE max=SIZE avg=SIZE
                // Read in values to memory
                x_i = x[j];
                y_i = y[j];
                y_err_i = y_err[j];

                // Calculate
                t_i = (x_i - S_x/S)/y_err_i;
                S_tt += t_i * t_i;
                
                b_partial += (t_i * y_i)/y_err_i;

                if (last[j] == 1){
                    b = b_partial / S_tt;
                    a = (S_y - S_x * b)/S;
                    
                    sigma_a2_partial = 1 + (S_x * S_x)/(S*S_tt);
                    sigma_a2 = sigma_a2_partial/S;
                    sigma_b2 = 1/S_tt;
                }
            }
            // Goodness of fit
            chi = 0;
            chi_2 = 0;
            chi2_ndof = 0;
            for (int k=event_start; k <= event_end; k++){
                #pragma HLS LOOP_TRIPCOUNT min=SIZE max=SIZE avg=SIZE
                chi = (y[k] - a - b*x[k])/y_err[k];
                chi_2 += chi * chi;
                if (last[k] == 1){
                    int events = event_end - event_start + 1;
                    event_start = event_end + 1;
                    std::cout << events << std::endl;
                    chi2_ndof = chi_2 / (events - 2);
                    



                    std::cout << "Event fit" << std::endl;
                    std::cout << "Slope: " << b << " +/- " << sqrt(sigma_b2) << std::endl;
                    std::cout << "Y intercept: " << a << " +/- " << sqrt(sigma_a2) << std::endl;
                    std::cout << "Chi_2/ndof: " << chi2_ndof << std::endl;
                    std::cout << std::endl;
                    
                    // Reset sums to 0
                    S_x = 0;
                    S_y = 0;
                    S = 0;
                }
            }
            
        }
    }
}
}

      

