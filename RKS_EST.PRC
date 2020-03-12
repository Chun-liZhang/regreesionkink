/*
RKS_EST.PRC

written by:

Lixiong Yang
School of Management, LZU

This procedure computes estimates and confidence
intervals for threshold models with state-dependent threshold.  
The procedure takes the form

the step could be reset following your data.

{q1hat,ups1hat,ups2hat}= rks_est(dat,names,yi,xi,qi,qi1,qi2,h,kk)

The inputs are:
  dat     = data matrix (nxk)
  names   = variable names (kx1), corresponding to dat matrix
  yi      = index of dependent (y) variable, e.g.: yi = 1
  xi      = indexes of independent (x) variables, e.g.: xi = 2|3
  qi      = index of threshold (q) variable, e.g.: qi = 4;
  h       = heteroskedasticity indicator
            Set h=0 to impose homoskedasticity assumption
            Set h=1 to use White-correction for heteroskedasticity
  kk      = frequency

Output:
  q1hat    = r0
  ups1hat=LS estimate of coefficent of sin
  ups2hat=LS estimate of coefficent of cos

time-varying threshold: r0+ups1hat*sin+ups2hat*cos

The remaining outputs are printed to the screen.

Notes:  
  (1)  Do not include a constant in the independent variables;
       the program automatically adds an intercept to the regression.

  (2)  There are three other control parameters, governing the choice
       of confidence level, the nonparametric method used to compute
       the nuisance parameter in the event of heteroskedastic, and 
       whether to print the graph of the likelihood.  These controls
       are listed at the beginning of the procedure code.

  (3)  The program includes code to produce a graph.
       This requires that the pgraph library is active.
       Thus, the command "library pgraph" must have been issued
       previously in the GAUSS session.
       if this is a problem, delete the lines between
          "if _graph==1;"  and  "endif;"



Example:
if the nxk matrix "dat" contains the dependent variable in the first
column, the independent variables in the second through tenth columns,
and the threshold variable in the fifth.  if the error is homoskedastic:

    xi = 2|3|4|5|6|7|8|9|10;
    qhat = thr_est(dat,1,xi,5,0);

while if the error is (possibly)heteroskedatic, 
replace the second line with

   {q1hat,ups1hat,ups2hat}= tv_thr_est(dat,1,xi,5,1,kk);


*******************************************************************/

proc (5)=rks_est(dat,names,yi,xi,qi,qi1,qi2,zi,h,kk);

local n,qs,q,y,x,k,yname,xname,qname,mi,beta,e,ee,sig,xe,se,vy,r_2,
qn,sn,irb,mm,sume,ci,r,irf,ir,mmi,rmin,smin,qhat,sighat,xx,
i1,i2,x1,y1,x2,y2,mi1,mi2,beta1,beta2,e1,e2,ej,n1,n2,ee1,ee2,sig1,sig2,
sig_jt,se1,se2,vy1,vy2,r2_1,r2_2,r2_joint,lr,r1,r2,qx,qh,m1,m2,g1,g2,
sigq,hband,u,u2,f,df,eps,sige,kh,eta2,c1,c2,lr1,lr2,qhat1,qhat2,z,
beta1l,beta1u,beta2l,beta2u,b1,b2,ser1,ser2,ff,p,tit,vtit,clr,plr,irr,upss1,upss2,rr0hat,rr1hat,rr2hat,
_conf1,_conf2,_nonpar,_graph,qqt,qqs,r0hat,r1hat,r2hat,mssr,g0hat,g1hat,g2hat,numbb,
ups1,rg1,ups1hat,ups2,rg2,ups2hat,nns,nreps,trd,q0hat,ssrh,rrmin,ssmin,q1hat,jj1,jj2,start1,step1,rgtv,
irff,x11,y11,x22,y22,mi11,mi22,beta11,beta22,e11,e22;

@ Control Parameters, can be modified if desired  @

_conf1 = .95;  @ Confidence Level for Confidence Regions  @
_conf2 = .8;   @ Confidence Level for first step of two-step
                 Confidence Regions for regression parameters @
_nonpar = 2;   @ Indicator for non-parametric method used to
                 estimate nuisance scale in the presence of
                 heteroskedasticity (only relevant if h=1).
                 Set _nonpar=1  to estimate regressions using
                 a quadratic.
                 Set _nonpar=2  to estimate regressions using
                 an Epanechnikov kernel with automatic bandwidth. @

_graph  = 1;   @ Set _graph=1 for the program to produce the graph
                 of the concentrated likelihood in gamma.
                 Set _graph=0 to not view the graph.             @

         @frequency of fourier terms@




if (h .ne 0)*(h .ne 1);
"You have entered h = " h;
"This number must be either 0 (homoskedastic case) or 1 (heteoskedastic)";
"The program will either crash or produce invalid results";
endif;
if (_nonpar .ne 1)*(_nonpar .ne 2)*(h==1);
"You have entered _nonpar = " _nonpar;
"This number should be either 1 (quadratic regression)"
"   or 2 (kernel regression)";
"The program will employ the quadratic regression method";"";"";
endif;

n = rows(dat);
q = dat[.,qi];
qs = sortind(q);/*??q¡ä¨®D?¦Ì?¡ä¨®??D¨°¦Ì???¡À¨º*/
q = q[qs];
y = dat[qs,yi];/*¨®?q¦Ì??3D¨°¨°???¦Ì???¨¢D?¨´¡À?*/
x = ones(n,1)~dat[qs,xi];
xx=dat[qs,xi];
z=dat[qs,zi];
k = cols(x);
yname = names[yi];
qname = names[qi];
xname = "Constant"|names[xi]|names[zi];/*=====*/

@generating sin and cos@
   trd=seqa(1,1,n); /*trend, 1 to n*/
   upss1=sin(2*pi*kk*trd/n);
   upss2=cos(2*pi*kk*trd/n);
@generating sin and cos@

ups1=dat[qs,qi1];
ups2=dat[qs,qi2];



x=x~z;/*====*/
mi = invpd(moment(x,0));/*moment(x,0)=x'x*/
beta = mi*(x'y);/*inv(x'x)x'y*/
e=y-x*beta;
ee=e'e;
sig = ee/(n-k);
xe = x.*e;
if h==0;
  se = sqrt(diag(mi)*sig);
else;
  se = sqrt(diag(mi*moment(xe,0)*mi));
endif;
vy = sumc((y - meanc(y)).^2);/*SST*/
r_2 = 1-ee/vy;


@computing min ssr@

start1=-4; /*-1*maxc(q)*/
step1=0.1;       @sin?¦Ì¨ºy¦Ì?¨¨??¦Ì?????e¨º??¦Ìstart1¡ê?2?3¡è,nns??¨ºy?Y@
nns=floor(-2*start1/step1);@number of rg1 divided@

rg1=seqa(-4,step1,nns);/*sin?¦Ì¨ºy¦Ì?¨¨??¦Ì?????e¨º??¦Ì-1¡ê?step=0.1,nns??¨ºy?Y*/
rg2=seqa(-6,step1,nns);




qqs=q;
qs = unique(qqs,1);/*Sorts (small to large) and removes duplicate elements from a vector*/
qs=seqa(10,1,70);
qn = rows(qs);
sn = zeros(qn,1);




irb = zeros(n,1);
mm=zeros(k,k);
sume = zeros(k,1);
ci = 0;


q0hat=zeros(nns,nns);
@fix r0, varying r1 r2@
r1hat=zeros(qn,1);r2hat=zeros(qn,1);ssrh=zeros(nns,nns);mssr=zeros(qn,1);

r=1; do while r<=qn;

  for j1(1,nns,1);
  for j2(1,nns,1);

    rgtv=rg1[j1]*ups1+rg2[j2]*ups2;  

  irf = (q .<= qs[r]+rgtv);/*?¨²¡Á?¨¬??t=1*/

  irr = irf - irb;/*¡Áa?¡¥?a???¨®?¨´¡À?*/
 
ci=sumc(irf);



irff = 1-irf;
x11 = (xx-qs[r]-rgtv).*irf;

x22 = (xx-qs[r]-rgtv).*irff;
y11=x11~x22~ones(n,1)~z;



  if (ci .> n/10)*(ci .< (n-n/10));/*?????¨²???D¦Ì?2?¨ºyD??¨¹¨º?¡Àe*/
    
  
 @SSR computed by new approach@

mi11 = invpd(moment(y11,0));

beta11 = mi11*(y11'y);

e11 = y - y11*beta11;

    @SSR computed by new approach @
    ssrh[j1,j2] = e11'e11;
   

  else;
   ssrh[j1,j2]=ee;
  endif;



 endfor;
endfor;

@mod: find the row and column of min@
rrmin=minindc(vecr(ssrh));/*the index of the smallest element in each column of sn*/
if fmod(rrmin,nns)==0;
jj1=rrmin/nns;
jj2=nns;
else;
jj2=fmod(rrmin,nns);/*colums,index of ups2*/
jj1=floor(rrmin/nns)+1;/*row,index of ups1*/
endif;
r1hat[r]=rg1[jj1];

r2hat[r]=rg2[jj2];
mssr[r]=ssrh[jj1,jj2];


  

r=r+1;endo;






@===@


@time-vary threshold can divided sample as the same way for diffrent threshold, which makes min ssr unsigle@
rrmin=minc(vecr(mssr));/*the index of the smallest element in each column of sn*/

numbb=0;
rr0hat=0;rr1hat=0;rr2hat=0;

for i(1,qn,1);
if (mssr[i]==rrmin);
rr0hat=rr0hat+qs[i];
rr1hat=rr1hat+r1hat[i];
rr2hat=rr2hat+r2hat[i];
numbb=numbb+1;

else;
numbb=numbb;
endif;

endfor;






rmin=minindc(mssr);/*the index of the smallest element in each column of sn*/

ssmin=mssr[rmin];/*min ssr*/
q1hat=qs[rmin];/*est of threshold*/
ups1hat=rr1hat/numbb;
ups2hat=rr2hat/numbb;
 


sighat=ssmin/n;

qqt=q1hat+ups1hat*ups1+ups2hat*ups2;



i1=(q.<=qqt);
i2=1-i1;

x1=(xx-qqt).*i1;

x2=(xx-qqt).*i2;
y1=x1~x2~ones(n,1)~z;
mi1=invpd(moment(y1,0));

beta1=mi1*(y1'y);


e1=y-y1*beta1;

n1=sumc(i1);

sig1=e1'e1/n;






"qs~r1hat~r2hat~mssr=" qs~r1hat~r2hat~mssr~(mssr-ssmin)*n/ssmin;

@trd=seqa(1,1,n)+1791;
library pgraph;
   graphset;
xy(qs,(mssr-ssmin)*n/ssmin);@

@xy(trd,qqt);@

@output result@

ff = "#*.*lG"~20~8;
"Global OLS Estimation, Without Threshold";
print ("Dependent Variable:     " $+ yname);
if h==1;
  "Heteroskedasticity Correction Used";
elseif h==0;
  "OLS Standard Errors Reported";
endif;
"";
p = printfm("Variable"~"Estimate"~"St Error",0~0~0,ff);"";
"       __________________________________________";
p = printfm(xname~beta~se,0~1~1,ff);"";
"";
"Observations                      " n;
"Degrees of Freedom                " (n-k);
"Sum of Squared Errors             " ee;
"SST="            vy;
"Residual Variance                 " sig;
"R-squared                         " r_2;
"Heteroskedasticity Test (P-Value) no";
"";"";
"__________________________________________________";



"Threshold model with time-varying value:";

"Frequency of fouriers kk=" kk;
"Min SSR of time-varying model=" ssmin;
"R2=SSR/SST" 1-ssmin/vy;


"Threshold Estimation";

"Threshold Estimate constant       " q1hat;
"Threshold q1" ups1hat;
"threshold q2" ups2hat;
"==========";

"parameters";


"Observations (x<=r_t)                     " n1;

"beta1-threshold=" beta1;

"sig1=" sig1;








retp(q1hat,ups1hat,ups2hat,beta1,ssmin);
endp;
