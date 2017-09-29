
# include <iostream>
# include <armadillo>

void fill_array(arma::mat& A, int n){
    double rho_0 = 0.0;
    double rho_n = 4.0;
    arma::vec rho(n+2);
    rho(0) = rho_0;
    rho(n+1) = rho_n;

    double h_step = (rho_n - rho_0)/(n+1);
    double hh =h_step*h_step;

    //std::cout << "h_step= " << h_step << std::endl;
    //std::cout << "hh = " << hh << std::endl;

    for (int i=1; i<n+1; i++){
        rho(i) = rho_0 + i*h_step;
    }

    //std::cout << "rho" << rho << std::endl;

    arma::vec diag_el(n+1);
    for (int i=0; i<n+1; i++){
        diag_el(i)= (2.0/hh) + (rho(i)*rho(i));
    }

    double off_const = -1.0/hh;

    //I am not including rho_0
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j){A(i,j)=diag_el(i+1);}
            if (fabs(i-j) == 1){A(i,j)=off_const;}

        }
    }
    //diag_el.print("Diag element");
}

void fill_array_interactive(arma::mat& A, int n){
    double rho_0 = 0.0;
    //Rho max scales with the frequency omega
    double rho_n = 7.0;
    arma::vec rho(n+2);
    rho(0) = rho_0;
    rho(n+1) = rho_n;

    double h_step = (rho_n - rho_0)/(n+1);

    std::cout << "h_step= " << h_step << std::endl;

    double hh =h_step*h_step;

    std::cout << "hh = " << hh << std::endl;

    for (int i=1; i<n+1; i++){
        rho(i) = rho_0 + i*h_step;
    }
    //std::cout << "rho= " << rho << std::endl;

    //Defining the frequency, 0.01, 0.5, 1, 5
    double omega = 0.500;
    double omega_squared = omega*omega;


    arma::vec diag_el(n);

    for (int i=0; i<n; i++){
        diag_el(i)= (2.0/hh) + (rho(i+1)*rho(i+1))*omega_squared + (1.0/rho(i+1));
    }

    //diag_el.print("Diag element = ");

    double off_const = -1.0/hh;

    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            if (i==j){A(i,j)=diag_el(i);}
            if (fabs(i-j) == 1){A(i,j)=off_const;}

        }
    }

    //A.print("A= ");

}

