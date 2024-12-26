#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lambda(NumericVector baselambda, NumericVector explp) {
  int n=baselambda.size();
  NumericVector result(n);
  result=result+baselambda(0);
  result[Range(1,n-1)]=baselambda[Range(1,n-1)]-baselambda[Range(0,n-2)];
  return result*explp;
}

// [[Rcpp::export]]
NumericVector S1(NumericVector baselambda12, NumericVector explp12,
                 NumericVector baselambda13,NumericVector explp13) {
  NumericVector lambda12=lambda(baselambda12, explp12);
  NumericVector lambda13=lambda(baselambda13, explp13);
  NumericVector temp=cumsum(log(1-pmin(lambda12+lambda13,0.99999)));
  return exp(temp);
}

// [[Rcpp::export]]
NumericVector S2(NumericVector baselambda23, NumericVector explp23) {
  NumericVector lambda23=pmin(lambda(baselambda23, explp23),0.99999);
  NumericVector temp=cumsum(log(1-lambda23));
  return exp(temp);
}

// [[Rcpp::export]]
NumericVector P11(NumericVector baselambda12, NumericVector explp12,
                  NumericVector baselambda13,NumericVector explp13) {
  return S1(baselambda12, explp12, baselambda13, explp13);
}

// [[Rcpp::export]]
NumericVector P122(NumericVector baselambda12,NumericVector explp12,
                   NumericVector baselambda13, NumericVector explp13,
                   NumericVector baselambda23, NumericVector explp23) {
  NumericVector lambda12=lambda(baselambda12, explp12);
  NumericVector survs=S1(baselambda12, explp12, baselambda13, explp13);
  int n=survs.size();
  survs.push_front(1);
  survs.erase(n);
  NumericVector survs2=S2(baselambda23,explp23);
  NumericVector result(n);
  for(int i=0; i < n; i++) {
    NumericVector temp=lambda12[Range(0,i)]*survs[Range(0,i)]*survs2(i)/survs2[Range(0,i)];
    result(i)=sum(temp);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector P123(NumericVector baselambda12,NumericVector explp12,
                   NumericVector baselambda13, NumericVector explp13,
                   NumericVector baselambda23, NumericVector explp23) {
  NumericVector lambda12=lambda(baselambda12, explp12);
  NumericVector survs=S1(baselambda12, explp12, baselambda13, explp13);
  int n=survs.size();
  survs.push_front(1);
  survs.erase(n);
  NumericVector survs2=S2(baselambda23,explp23);
  NumericVector result(n);
  for(int i=0; i < n; i++) {
    NumericVector temp=lambda12[Range(0,i)]*survs[Range(0,i)]*(1-survs2(i)/survs2[Range(0,i)]);
    result(i)=sum(temp);
  }
  return result;
}

// [[Rcpp::export]]
NumericVector P113(NumericVector baselambda12,NumericVector explp12,
                   NumericVector baselambda13,NumericVector explp13) {
  NumericVector lambda13=pmin(lambda(baselambda13, explp13),0.99999);
  NumericVector survs=S1(baselambda12, explp12, baselambda13, explp13);
  int n=survs.size();
  survs.push_front(1);
  survs.erase(n);
  return cumsum(lambda13*survs);
}

// [[Rcpp::export]]
List P(NumericVector times,
       NumericVector basehaz12,NumericMatrix explp12,
       NumericVector basehaz13,NumericMatrix explp13,
       NumericVector basehaz23,NumericMatrix explp23) {
  
  int n=explp12.nrow();
  int n_t=times.size();
  List result=List::create();
  for(int i=0; i < n; i++) {
    NumericMatrix p_i(3,n_t);
    NumericVector explp12_i(n_t);
    NumericVector explp13_i(n_t);
    NumericVector explp23_i(n_t);
    explp12_i=explp12_i+explp12(i,_);
    explp13_i=explp13_i+explp13(i,_);
    explp23_i=explp23_i+explp23(i,_);
    p_i(0,_)=P11(basehaz12,explp12_i,basehaz13,explp13_i);
    p_i(1,_)=P122(basehaz12,explp12_i,basehaz13,explp13_i,basehaz23,explp23_i);
    p_i(2,_)=P113(basehaz12,explp12_i,basehaz13,explp13_i)+P123(basehaz12,explp12_i,basehaz13,explp13_i,basehaz23,explp23_i);
    result.push_back(p_i);
  }
  return result;
}

// [[Rcpp::export]]
double D(NumericVector x_i,NumericVector x_j,NumericVector x_k,NumericVector p_i,NumericVector p_j, NumericVector p_k) {
  return sqrt(sum(pow(x_i-p_i,2)))+sqrt(sum(pow(x_j-p_j,2)))+sqrt(sum(pow(x_k-p_k,2)));
}

// [[Rcpp::export]]
double concordance_result(NumericVector p_i, NumericVector p_j,  NumericVector p_k) {
  
  NumericMatrix obs(3,3);
  obs(0,0)=1;
  obs(1,1)=1;
  obs(2,2)=1;
  
  NumericMatrix perms(6,3);
  NumericVector col1={0,0,1,1,2,2};
  NumericVector col2={1,2,0,2,0,1};
  NumericVector col3={2,1,2,0,1,0};
  perms(_,0)=col1;
  perms(_,1)=col2;
  perms(_,2)=col3;
  
  NumericVector D_vec(6);
  for(int i=0; i < 6; i++) {
    D_vec(i)=D(obs(_,perms(i,0)),obs(_,perms(i,1)),obs(_,perms(i,2)),p_i,p_j,p_k);
  }
  return which_min(D_vec)==0;
}

// [[Rcpp::export]]
List VUS_ID_t(double t, int t_index, NumericVector ill_times, NumericVector death_times, NumericVector censoring_times,
              List transition_probs) {
  
  NumericVector ind(death_times.size());
  for(int i=0; i < ind.size(); i++) {
    ind(i)=i;
  }
  
  NumericVector ind_initial;
  NumericVector ind_ill;
  NumericVector ind_death;
  
  double n_conc=0;
  double n_comp=0;
  
  ind_initial=ind[ill_times > t & death_times > t & censoring_times > t];
  ind_ill=ind[ill_times < t & death_times > t & censoring_times > t];
  ind_death=ind[death_times==t & censoring_times > t];
  int n_initial=ind_initial.size();
  int n_ill=ind_ill.size();
  int n_death=ind_death.size();
  
  if(n_initial > 0 & n_ill > 0 & n_death >0) {
    for(int i=0; i < n_initial; i++) {
      for(int j=0; j < n_ill; j++) {
        for(int k=0; k < n_death; k++) {
          NumericMatrix transition_probs_i=transition_probs[ind_initial(i)];
          NumericMatrix transition_probs_j=transition_probs[ind_ill(j)];
          NumericMatrix transition_probs_k=transition_probs[ind_death(k)];
          n_conc=n_conc+concordance_result(transition_probs_i(_,t_index),
                                           transition_probs_j(_,t_index),
                                           transition_probs_k(_,t_index));
        }
      }
    }
  }
  
  n_comp=n_comp+n_initial*n_ill*n_death;
  
  return List::create(n_conc/n_comp,n_conc,n_comp);
}

// [[Rcpp::export]]
List VUS_CD_t(double t, int t_index, NumericVector ill_times, NumericVector death_times, NumericVector censoring_times,
              NumericVector G,List transition_probs) {
  
  NumericVector ind(death_times.size());
  for(int i=0; i < ind.size(); i++) {
    ind(i)=i;
  }
  
  NumericVector ind_initial;
  NumericVector ind_ill;
  NumericVector ind_death;
  
  double n_conc=0;
  double n_comp=0;
  
  ind_initial=ind[ill_times > t & death_times > t & censoring_times > t];
  ind_ill=ind[ill_times < t & death_times > t  & censoring_times > t];
  ind_death=ind[death_times < t  & censoring_times > t];
  int n_initial=ind_initial.size();
  int n_ill=ind_ill.size();
  int n_death=ind_death.size();
  
  if(n_initial > 0 & n_ill > 0 & n_death >0) {
    for(int i=0; i < n_initial; i++) {
      for(int j=0; j < n_ill; j++) {
        for(int k=0; k < n_death; k++) {
          NumericMatrix transition_probs_i=transition_probs[ind_initial(i)];
          NumericMatrix transition_probs_j=transition_probs[ind_ill(j)];
          NumericMatrix transition_probs_k=transition_probs[ind_death(k)];
          n_conc=n_conc+(1/G(ind_death(k)))*concordance_result(transition_probs_i(_,t_index),
                                            transition_probs_j(_,t_index),
                                            transition_probs_k(_,t_index));
          
          n_comp=n_comp+(1/G(ind_death(k)));
        }
      }
    }
  }
  
  return List::create(n_conc/n_comp,n_conc,n_comp);
}

// [[Rcpp::export]]
List VUS_ID(NumericVector times,
         NumericVector ill_times, NumericVector death_times, 
         NumericVector censoring_times,
         NumericVector basehaz12, NumericMatrix explp12,
         NumericVector basehaz13, NumericMatrix explp13,
         NumericVector basehaz23, NumericMatrix explp23) {
  
  List transition_probs=P(times,basehaz12,explp12,basehaz13,explp13,basehaz23,explp23);
  
  int n=times.size();
  NumericVector result(n);
  NumericVector conc(n);
  NumericVector comp(n);
  for(int i=0; i < n; i++) {
    /*Rcout << i << std::endl;*/
    List vus_t=VUS_ID_t(times(i),i,ill_times,death_times,censoring_times,transition_probs);
    result(i)=vus_t(0);
    conc(i)=vus_t(1);
    comp(i)=vus_t(2);
  }
  NumericVector result_clean=result[!is_nan(result)];
  NumericVector conc_clean=conc[!is_nan(result)];
  NumericVector comp_clean=comp[!is_nan(result)];
  NumericVector times_clean=times[!is_nan(result)];
  return List::create(result_clean,conc_clean,comp_clean,times_clean);
}

// [[Rcpp::export]]
List VUS_CD(NumericVector times,
            NumericVector ill_times, NumericVector death_times,
            NumericVector censoring_times,NumericVector G,
            NumericVector basehaz12, NumericMatrix explp12,
            NumericVector basehaz13, NumericMatrix explp13,
            NumericVector basehaz23, NumericMatrix explp23) {
  
  List transition_probs=P(times,basehaz12,explp12,basehaz13,explp13,basehaz23,explp23);
  
  int n=times.size();
  NumericVector result(n);
  NumericVector conc(n);
  NumericVector comp(n);
  for(int i=0; i < n; i++) {
    Rcout << i << std::endl;
    List vus_t=VUS_CD_t(times(i),i,ill_times,death_times,censoring_times,G,transition_probs);
    result(i)=vus_t(0);
    conc(i)=vus_t(1);
    comp(i)=vus_t(2);
  }
  NumericVector result_clean=result[!is_nan(result)];
  NumericVector conc_clean=conc[!is_nan(result)];
  NumericVector comp_clean=comp[!is_nan(result)];
  NumericVector times_clean=times[!is_nan(result)];
  return List::create(result_clean,conc_clean,comp_clean,times_clean);
}




