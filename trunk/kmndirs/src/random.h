/* header for random.c */
double   runi(void);
double runir(double a,double b);
int runii(int na,int nb);
double rnor(double mu,double sd);
int rpois(double mu);
double rexp(double lambda);
double rcauchy(double loc,double scale);
double rncauchy(double loc,double scale);
double rgamma(double alpha);
double rbeta(double alpha,double beta);
double rlogistic(double loc,double scale);
double rlnorm(double norm,double norsd);
int rbin(int n,double p);
double rweibull(double gamma);
double rchisq(double t);
double rf(double t,double u);
double rstudent(double t);
void  setseed(unsigned int s);
long  genseed(void);


