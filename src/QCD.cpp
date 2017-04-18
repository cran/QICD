//' @export for cpp
#include <iostream>
#include <math.h>
#include <R.h>
using namespace std;
//extern "C" __declspec(dllexport)
/*  
function beta_iteration for updating the beta. 
x is the predictor matrix, y is the reponse viable. 
beta1 is the coefficients for updating.
j is the subscript of the excluded column of x
tao is the quantile parameter.
return the vector prebeta for caculating the weighted median 
*/
//extern "C"{
//function for iteration to get the optimal solution. 
// y is the responses, x is the predictors, beta1 is the parameters, nyrow is the row length 
// of y, nxcol is the column length of x, tao is the quantile parameter.
extern "C"{
void F77_NAME(xssort)(double*, double*, int*, int*);

// C WRAPPER FUNCTION FOR FORTRAN SUBROUTINE
void quicksort(double *x, double *w, int *n){
	int kflag = 2;
//  CALLING FORTRAN	FUNCTION FROM C
	F77_CALL(xssort)(x, w, n, &kflag);
}
void QCD(double *y, double *x, double *beta0, double *beta1,double *pre_value, int *nyrow, 
					int *nxcol, double *tao,double *lambda, double *a, int *index,int *index1,double *thresh,int *maxin)
//y is the response variable, nx1 vector; x is the predictor, nx(p+1) matrix; beta0 is the initial value of coefficients, (p+1)x1 vector ;
//beta1 will be the final value of coefficients, (p+1)x1 vector; pre_value is the vector to calculate the weighted median,nx1 vector,x%*%beta;
//nyrow is the length of y, which is n; nxcol is the number of columns of x, which is p+1; tao is the quantile number; lambda is the tuning 
//parameter; a is another parameter for SCAD and MCP penalties; index is the indicator for penalty types, 0 for SCAD, 1 for MCP and 2 for LASSO.
//index1 is the intercept indicator, 1 for including intercept, 0 for without intercept
//thresh is the convergence threshold
//maxin is maximal iteration numbers
{	
	
	double *pre_value_final=new double[*nyrow+1];
	// pre_value_final is the vector to save the weighted median vector
	double *weight=new double[*nyrow+1];
	// weight
	int length_weight=0;
	//weight length
	int count_xj=0;
	//length of nonezero x elements
	int countj;
	// index for the count row count5 column entry of x
	int count,count1,count3,count4,count5,end;
	double temp3,temp_sum,temp1,distance=0;
	count4=0;
	count1=0;
	end=(*nxcol)-1;
	while(count1<(*maxin))
	//number of outside iterations
	{	
		count5=count4%(*nxcol);
		//count5 column
		temp_sum=0;
		count_xj=0;
		//index of weight
		temp1=0;
		//count 5 is fix suffix of beta.
		for (count=0,countj=(*nyrow)*count5;count<*nyrow;count++,countj++)
		// count is the number for the yrow
		{	
			if (x[countj]!=0)
			//x[countj] corresponds to beta1x[count5]
			{
				//int previous=((*nxcol+count5-1)%(*nxcol));
				int previous=0;
			  if(count5==0)
			    previous=end;
			  else
			    previous=count5-1;
				pre_value[count]+=distance*x[count+(*nyrow)*previous];		
				//pre_value update
				pre_value_final[count_xj]=y[count]-pre_value[count];
				//weight terms.
				pre_value_final[count_xj]>0 ? weight[count_xj]=fabs(x[countj])*(*tao)/double(*nyrow): weight[count_xj]=
					fabs(x[countj])*(1-*tao)/double(*nyrow);
				temp_sum+=weight[count_xj];
				if ((*tao)>=0.5)
					pre_value_final[count_xj]=-(pre_value_final[count_xj]+x[countj]*beta1[count5])/double(x[countj]);
				else
					pre_value_final[count_xj]=(pre_value_final[count_xj]+x[countj]*beta1[count5])/double(x[countj]);
				count_xj++;
			}
		}
			if ((*index1)!=1||count5!=(*nxcol-1))
			{
				pre_value_final[count_xj]=0;

				/////////////////////////////////////////////////SCAD
				if ((*index)==0)
				{
      				if (fabs(beta0[count5])<(*lambda))
						weight[count_xj]=*lambda;
					else if (fabs(beta0[count5])<(*a)*(*lambda))
						weight[count_xj]=((*a)*(*lambda)-fabs(beta0[count5]))/(double(*a-1));
					else
						weight[count_xj]=0;
				}
				///////////////////////////////////////////////
				else if ((*index)==1)
				{
				//////////////////////////////////////////MCP
					if (fabs(beta0[count5])<(*a)*(*lambda))
						weight[count_xj]=(*lambda)-fabs(beta0[count5])/double(*a);
					else
						weight[count_xj]=0;
				///////////////////////////////////////////
				}
				//////////////////////////////////////////
				else if ((*index) == 2)
				{
					//////////////////////////////////////////LASSO
					weight[count_xj] = (*lambda);
					///////////////////////////////////////////
				}
				temp_sum+=weight[count_xj];
				length_weight=count_xj+1;
				quicksort(pre_value_final,weight,&length_weight);
				for (count3=0;count3<=count_xj;count3++)
				{
					temp1+=weight[count3];
					if (temp1> temp_sum/2)
					{
						break;
					}
				}
			}
			else
			{
				length_weight=count_xj;
				quicksort(pre_value_final,weight,&length_weight);
				for (count3=0;count3<count_xj;count3++)
				{
					temp1+=weight[count3];
					if (temp1> temp_sum/2)
					{
						break;
					}
				}
			}
			if ((*tao)>=0.5)
				temp3=-pre_value_final[count3];
			else
				temp3=pre_value_final[count3];
			//change from current beta1 to new beta1
			//distance=temp3-beta1[count5];
			//fabs(temp3)>(*thresh) ? beta1[count5]=temp3+(*thresh) : beta1[count5]=0;
			if(fabs(temp3)>(*thresh))
			{
			  distance=temp3-beta1[count5]+(*thresh);
			  beta1[count5]+=distance;
			}
			else
			{
			  distance=-beta1[count5];
			  beta1[count5]=0;
			}
			if ((count1==0||(fabs(distance)>(*thresh))||distance==0)&&count4<(*nxcol))
			{	
				count4++;
			}
			else 
			{	
				count4=0;
				count1++;
				end=count5;
				//fabs(distance)<=(*thresh) ? distance=0 : distance=distance;
				if(fabs(distance)<=(*thresh))
				{
				  beta1[count5]-=distance;
				  distance=0;
				}
			}
		
	}
		
		//*tao=distance;	
		delete [] pre_value_final;
		delete [] weight;
}
}
