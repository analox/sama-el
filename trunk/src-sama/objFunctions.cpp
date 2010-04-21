#define MAX_DIM 50
double rotation[MAX_DIM][MAX_DIM];
double translation[MAX_DIM];
unsigned dimension=30;

// Weierstrass global variable
double weierstrass_pa[100];
double weierstrass_pb[100];
double weierstrass_a = 0.5;
double weierstrass_b = 3;
double weierstrass_k = 20;
double weierstrass_sum2 = 0;
double weierstrass_pi = acos(-1.0);
// END: Weierstrass global variable

// Schwefel206
double schwefel206_A[100][100];
double schwefel206_B[100];
// END Schwefel206

// Schwefel213
double schwefel213_A[100][100];
double schwefel213_B[100][100];
double schwefel213_alpha[100];
// END Schwefel213
void initTransform()
{
	memset(rotation, 0, sizeof(rotation));
	memset(translation, 0, sizeof(translation));

	for (int i = 0; i < MAX_DIM; i++)
		rotation[i][i] = 1;

        // Weierstrass global variable
        for( int j = 0; j <= weierstrass_k; j++ ) {
                weierstrass_pa[j] = pow( weierstrass_a, j );
                weierstrass_pb[j] = pow( weierstrass_b, j );
                weierstrass_sum2 += weierstrass_pa[j] * cos(weierstrass_pi*weierstrass_pb[j]);
    	}
	// END Weierstrass global variable

	// Schwefel206
	int i, j;
        ifstream f("input_data/schwefel_206_data.txt");

        double schwefel206_o[100];
        for(i=0; i<100; i++)
                f >> schwefel206_o[i];

        // After loading the data, we set o[1..D/4] = -100 o[3D/4..D] = 100
        for(i=0; i<ceil(1.0*dimension/4); i++) schwefel206_o[i] = -100;
        for(i=floor(3.0*dimension/4)-1; i<dimension; i++) schwefel206_o[i] = 100;


        for(i=0; i<100; i++)
        {               
                for(j=0; j<100; j++)
                {
                        f >> schwefel206_A[i][j];
                }
        }


        for(i=0; i<dimension; i++)
        {               
                schwefel206_B[i] = 0;
                for(j = 0; j < dimension; j++)
                {
                        schwefel206_B[i] += schwefel206_A[i][j]*schwefel206_o[j];
                }
        }
	f.close();
	// END Schwefel206

	// Schwefel213
        ifstream f1("input_data/schwefel_213_data.txt");

        for(i=0; i<100; i++)
        {               
                for(j=0; j<100; j++)
                {
                        f1 >> schwefel213_A[i][j];                   
                }
        }

        for(i=0; i<100; i++)
        {               
                for(j=0; j<100; j++)
                {
                        f1 >> schwefel213_B[i][j];                   
                }
        }

        for(i=0; i<100; i++)
        {
                        f1 >> schwefel213_alpha[i];
        }

        f1.close();
	// END Schwefel213
}

void loadRotation(char * filename, int ndim)
{
	ifstream fin(filename);
	cout << filename << endl;
	for (int i = 0; i < ndim; i++) {
		for (int j = 0; j < ndim; j++)
		{
			fin >> rotation[i][j];
			cout << rotation[i][j] << " ";
		}
		cout << endl;
	}
	fin.close();
	
}

void loadTranslation(char * filename, int ndim)
{
	ifstream fin(filename);
	cout << filename << endl;
	for (int i = 0; i < ndim; i++) {
		fin >> translation[i];
		cout << translation[i] << " ";
	}
	cout << endl;
	fin.close();
}

vector<double> transform( vector<double> & x )
{
	int i, j;
	double x_[1000];
	vector<double> res(x.size(), 0);

	int ndim = x.size();	
	for( i = 0; i < ndim; i++ ) {
		x_[ i ] = ( x[ i ] - translation[ i ] );
	}

	for( i = 0; i < ndim; i++ ) {		
		for( j = 0; j < ndim; j++ ) {
			res[ i ] += rotation[ j ][ i ] * x_[ j ];
		}
	}
	return res;
}

// Sphere
double sphere( vector<double>& x_ ) {

	//
	vector<double>x = transform(x_);
	//
	unsigned size = x.size();
	unsigned i;
	double result=0;
	
	for(i=0; i<size; i++) {
		result += x[i]*x[i];
	
	}
	
	return result;
	
}

// Ackley
double ackley( vector< double >& x_ )
{
	//
	vector<double>x = transform(x_);
	//
    	const double A = 20.;
    	const double B = 0.2;
    	const double C = znPi2C;

    	unsigned i, n;
    	double   a, b;

    	for( a = b = 0., i = 0, n = x.size( ); i < n; ++i ) {
       		a += x[ i ] * x[ i ];
		b += cos( C * x[ i ] );
	    }

    	return -A * exp( -B * sqrt( a / n ) ) - exp( b / n ) + A + znEC;
}

// Rastrigin
double rastrigin( vector < double >& x_ ) {
	//
	vector<double>x = transform(x_);
	//
        int varsize = x.size( );
        int w, i, j, a = 2;
        float c;
        double result, unNormVar;

        result = 0;
        for ( w = 0; w < varsize; w++ ) {
                result = result + ( ( ( x[ w ] ) * ( x[ w ] ) ) - 10.0 * cos( 2 * M_PI * ( x[ w ] ) ) );
        }
        result = ( result + 10.0 * varsize );
        return result;
}

// Griewank
double griewank(vector<double>& x_)
{
	//
	vector<double> x = transform(x_);
	//
        int varsize = x.size();

        double l_value, l_Sumobj, l_Productobj;
        int w;

        l_value = 0;
        l_Sumobj = 0;
        l_Productobj = 1;


        for (w = 0; w < varsize; w++) {

                l_Sumobj = l_Sumobj + ((x[w] * x[w])/4000);
                l_Productobj = l_Productobj * cos(x[w]/sqrt(w+1.0));
        }
        l_value = (l_Sumobj + 1 - l_Productobj) ;

        return l_value;
}

// Bump
double bump( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;

        // check constraints
        double sumx = 0, prodx = 1;
        for (i=0; i<nDim; i++)
        {
                sumx += x[i];
                prodx *= x[i];
        }

        if ((sumx >= 7.5 * nDim) || (prodx <= 0.75)) return 0; // penalty


        double sumc4=0, prodc2=1, sumixi2=0;
        double tmp;

        for( i = 0; i < nDim; i++ ) {
                tmp = cos(x[i]);
                sumc4 += tmp*tmp*tmp*tmp;;
                prodc2 *= tmp*tmp;
                sumixi2 += (i+1)*x[i]*x[i];
        }
        
        return fabs(sumc4 - 2*prodc2) / sqrt(sumixi2);
}


// Elipptic
double elliptic( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;

        double res = 0.0;
        for( i=0; i<nDim; i++)
        {
                res += pow(1.0e6, (double)i/(nDim-1)) * x[i] * x[i];
        }

        return res;
}

// FGriewankRosenbrock
double griewank_rosenbrock( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;  
        
        double res = 0.0;       
        double* _x = new double[nDim];

        for (i=0; i<nDim; i++) _x[i] = x[i] + 1;

    	for (i=0; i<nDim-1; i++)
   	{
        	double temp = 100.0*pow((_x[i]*_x[i]-_x[i+1]),2.0) + pow((_x[i]-1.0),2.0);
                	res += (temp*temp)/4000.0 - cos(temp) + 1.0;
    	}

        double temp = 100.0*pow((_x[nDim-1]*_x[nDim-1]-_x[0]),2.0) + pow((_x[i]-1.0),2.0);
    	res += (temp*temp)/4000.0 - cos(temp) + 1.0;

    	return res; 
}
// Rosenbrock
double rosenbrock( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;
        double res;
        double _x[1000];

        for( i=0; i<nDim; i++) _x[i] = x[i] + 1;        

        res = 0;
        
        for( i = 0; i < nDim - 1; i++ ) {
                res += ( 100 * pow( _x[ i] * _x[ i] - _x[ i + 1], 2 ) + pow( _x[ i] - 1, 2 ) );
        }

        return res;
}

double scaffer( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;

        double res = 0.0;

    	for (i=0; i<nDim-1; i++)
    	{
       	 	double temp1 = pow((sin(sqrt(pow(x[i],2.0)+pow(x[i+1],2.0)))),2.0);
        	double temp2 = 1.0 + 0.001*(pow(x[i],2.0)+pow(x[i+1],2.0));
        	res += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
   	}
    	double temp1 = pow((sin(sqrt(pow(x[nDim-1],2.0)+pow(x[0],2.0)))),2.0);
    	double temp2 = 1.0 + 0.001*(pow(x[nDim-1],2.0)+pow(x[0],2.0));

    	res += 0.5 + (temp1-0.5)/(pow(temp2,2.0));
    	return res; 
}

double schwefel102( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i, j;
                
        double res = 0.0;

        for( i=0; i<nDim; i++)
        {
                double tmp = 0;
                for( j=0; j<=i; j++)
                {
                        tmp += x[j];
                }
                res += tmp*tmp;
        }

        return res;
}

double step( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
        int i;
        double sum1;
        
        sum1 = 0;       

        for( i = 0; i < nDim; i++ ) {
                sum1 += floor(x[i]);
        }

        return sum1 + 6*nDim;
}

double weierstrass( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//
	
        double sum2 = weierstrass_sum2 * nDim;

        int i, j;
        double sum1;

	sum1 = 0;

        for( i = 0; i < nDim; i++ ) {
                for( j = 0; j <= weierstrass_k; j++ ) {
                        sum1 += (  weierstrass_pa[j] * cos( 2 * weierstrass_pi * weierstrass_pb[j] * (x[ i ] + 0.5) ) );
                }
        }

        return sum1 - sum2;
}

//--------------------------------
        
double schwefel206( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//

	int i, j;
        double res = INT_MIN;

        
        for( i=0; i< nDim; i++)
        {
                double tmp = 0.;
                for(j=0; j<nDim; j++)
                {
                        tmp+=schwefel206_A[i][j]*x[j];
                }               
                tmp = fabs(tmp-schwefel206_B[i]);
                if (res < tmp) res = tmp;
        }       

        return res;
}

//------------------------------------------------
double schwefel213( vector<double>& x_ )
{
	//
	vector<double> x = transform(x_);
	int nDim = x.size();
	//

	int i, j;
        double res = 0.0;
        
        for( i=0; i< nDim; i++)
        {
                double sum1 = 0.0;
        	double sum2 = 0.0;
        	for (j=0; j<nDim; j++)
	        {
        		sum1 += schwefel213_A[i][j]*sin(schwefel213_alpha[j]) + schwefel213_B[i][j]*cos(schwefel213_alpha[j]);
	        	sum2 += schwefel213_A[i][j]*sin(x[j]) + schwefel213_B[i][j]*cos(x[j]);
        	}
        	res += pow((sum1-sum2),2.0);
        }

        return res;
}


