library(Rcpp) # Attach R package "Rcpp"
#find the index of min
cppFunction('List findmin(NumericVector x,int l){
double min=x[0];
int k;
for(int i=1;i<l;++i){
if(min>x[i]) {min=x[i];k=i+1;}
}
List s=List::create(min,k);
return s;
}')

