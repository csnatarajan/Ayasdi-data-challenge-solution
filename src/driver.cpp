#include "all.hpp"
#include "dlib/svm.h"
typedef double  ST;
typedef int ordinal;
using namespace std;
using namespace dlib;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::TimeMonitor;
using Teuchos::Time;
typedef matrix<ST, 310, 1> vec_space;
typedef radial_basis_kernel<vec_space> ktype;
typedef decision_function<ktype> func;
using std::vector;

string TableHelp::line() const
{
        ostringstream oss;
        for(ordinal i = 0; i < pagewidth_; i++)
        {
                oss << "-";
        }
        return oss.str();
}

string TableHelp::bigline() const
{
        ostringstream oss;
        for(ordinal i = 0; i < pagewidth_; i++)
        {
                oss << "=";
        }
        return oss.str();
}

string TableHelp::starline() const
{
        ostringstream oss;
        for(ordinal i = 0; i < pagewidth_; i++)
        {
                oss << "*";
        }
        return oss.str();
}

template<class T>
void write_data(const MAT2D<T> & Data, const string name)
{
        fstream f;
        f.open(name.c_str(), ios::out | ios::binary);
	if (! f)
        {
                cout << "Could not open file " << name << endl;
                exit;
        }
	f.write(reinterpret_cast<char*>(&Data.data[0][0]), sizeof(T)*Data.sizex*Data.sizey);
        f.close();
}

template<class T>
void write_data_asc(const MAT2D<T> & Data, const string & name)
{
	fstream f; f.open(name.c_str(), ios::out);
	for(ordinal i = 0; i < Data.sizey; i++)
	{
		for(ordinal j = 0; j < Data.sizex; j++)
			f << Data.data[i][j] << " ";
		f << endl;
	}
	f.close();
}

template<class T>
void write_data1d_asc(const MAT1D<T> & Data, const string & name)
{
        fstream f; f.open(name.c_str(), ios::out);
        for(ordinal i = 0; i < Data.sizex; i++)
        {
		f << Data.data[i] << endl;
        }
        f.close();
}


template<class T>
void get_mean(const MAT2D<T> & Data, MAT2D<T> & Avgdata, MAT1D<T> & Sum)
{
	for(ordinal i = 0; i < Sum.sizex; i++)
		Sum.data[i] /= Data.sizey;
	for (ordinal i = 0; i < Data.sizey; i++)
		for(ordinal j = 0; j < Data.sizex; j++)
			Avgdata.data[i][j] = Data.data[i][j] - Sum.data[j];
}

template<class T>
void get_mean2d(const MAT2D<T> & Data, MAT2D<T> & Avgdata)
{
	MAT1D<ST> Sum(Data.sizex);std::fill(Sum.data, Sum.data+Sum.sizex, 0.);
	for (ordinal i = 0; i < Data.sizey; i++)
	        for(ordinal j = 0; j < Data.sizex; j++)
        	        Sum.data[j] += Data.data[i][j];
        for (ordinal i = 0; i < Data.sizey; i++)
                for(ordinal j = 0; j < Data.sizex; j++)
                        Avgdata.data[i][j] = Data.data[i][j] - (Sum.data[j]/Data.sizey);
}

 
template<class T>
void get_mean1d(const MAT1D<T> & Data, MAT1D<T> & Avgdata)
{               
	ST sum;
        for(ordinal i = 0; i < Data.sizex; i++)
                sum += Data.data[i];
	sum /= Data.sizex;
        for (ordinal i = 0; i < Data.sizex; i++)
		Avgdata.data[i] = Data.data[i] - sum;
} 

pair<ST, ST> minmax(const MAT2D <ST> & Data)
{
pair<ST, ST> values;
ordinal i, N; 
	N = Data.sizex*Data.sizey;
	if (N%2 == 0)
	{        
		if (Data.space[0] > Data.space[1])    
		{
			values.second = Data.space[0];
			values.first = Data.space[1];
		} 
		else
		{
			values.first = Data.space[0];
			values.second = Data.space[1];
		}
		i = 2;
	} 
	else
	{
		values.first = Data.space[0];
		values.second = Data.space[0];
		i = 1;  
	}
	while (i < N-1) 
	{         
		if (Data.space[i] > Data.space[i+1])         
		{
			if(Data.space[i] > values.second)       
				values.second = Data.space[i];
			if(Data.space[i+1] < values.first)         
				values.first = Data.space[i+1];
		}
		else        
		{
			if (Data.space[i+1] > values.second)  
				values.second = Data.space[i+1];
			if (Data.space[i] < values.first)         
				values.first = Data.space[i];
		}       
		i += 2; 
	}           
return values;
}   
 
template<class T>
void histogram2(const MAT2D<T> & Data, const ordinal & bins, const string & fname)
{
	ordinal min, max;
	pair<ST, ST> values;
	values = minmax(Data);
	MAT1D<ordinal> hist(bins); fill(hist.data, hist.data+hist.sizex, 0);
	ST temp = (1.*bins)/(values.second - values.first);
	ordinal half = (-1*values.first*temp)+1;
	fstream f; f.open(fname.c_str(), ios::out);
	for(ordinal i = 0; i < Data.sizey; i++)
		for(ordinal j = 0; j < Data.sizex; j++)
		{
			if (Data.data[i][j] < 0)
				hist.data[half-(ordinal)floor(fabs(Data.data[i][j]*temp))] += 1;
			else
				hist.data[half+(ordinal)ceil(fabs(Data.data[i][j]*temp))] += 1;
		}
        for (ordinal i = 0; i < 200; i++)
        	f << fixed << values.first+i/temp << "\t" << hist.data[i] << endl;
        f.close();
}

template<class T>
void histogram1(const MAT1D<T> & Data, const ordinal & bins)
{
        ordinal min, max; min = -2; max = 2;
        MAT1D<ordinal> hist(bins); fill(hist.data, hist.data+hist.sizex, 0);
        ST temp = (1.*bins)/(max - min);
        ordinal half = bins/2;
	inputstringstream name("histogram");
        fstream f; f.open(name.c_str(), ios::out);
	for(ordinal i = 0; i < Data.sizex; i++)
	{
		if (Data.data[i] < 0)
			hist.data[half-(ordinal)floor(fabs(Data.data[i]*temp))] += 1;
		else
			hist.data[half+(ordinal)ceil(fabs(Data.data[i]*temp))] += 1;
	}
        for (ordinal i = 0; i < 200; i++)
                f << fixed << min+i/temp << "\t" << hist.data[i] << endl;
        f.close();
}

void get_variance(const MAT2D<ST> & Data, const std::string & name)
{
	MAT1D<ST> variance(Data.sizex); std::fill(variance.data, variance.data+variance.sizex, 0.);
	for(ordinal i = 0; i < Data.sizey; i++)
		for(ordinal j = 0; j < Data.sizex; j++)
			variance.data[j] += Data.data[i][j]*Data.data[i][j];
	ofstream odatafile(name.c_str());
        for(ordinal i = 0; i < variance.sizex; i++)
        	odatafile << variance.data[i]/(Data.sizey-1) << endl;
        odatafile.close();	
}


void get_coeffv(const MAT2D<ST> & Data, const std::string & name)
{
        MAT1D<ST> variance(Data.sizex), sum(Data.sizex);
	std::fill(variance.data, variance.data+variance.sizex, 0.); std::fill(sum.data, sum.data+sum.sizex, 0.);
	
        for(ordinal i = 0; i < Data.sizey; i++)
                for(ordinal j = 0; j < Data.sizex; j++)
		{
                        variance.data[j] += Data.data[i][j]*Data.data[i][j];
			sum.data[j] += Data.data[i][j];
		}

        ofstream odatafile(name.c_str());
        for(ordinal i = 0; i < variance.sizex; i++)
                odatafile << sqrt(variance.data[i]/(Data.sizey-1))/(sum.data[i]/Data.sizey) << endl;
        odatafile.close();
}

void get_cv_and_eig(const MAT2D<ST> & Data, MAT1D<ST> & evalues, MAT2D<ST> & evectors)
{
	MAT2D<ST> covarray(Data.sizex, Data.sizex);
	MAT1D<ordinal> suppz(covarray.sizex);
	ST vl, vu, sum = 0;
	ordinal N, il, iu, numvalues; il = covarray.sizex-1; iu = covarray.sizex;

	std::fill(covarray.space, covarray.space+(covarray.sizex*covarray.sizey), 0.);
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, covarray.sizey, covarray.sizex, Data.sizey,1, Data.space, Data.sizex, Data.space, Data.sizex, 0., covarray.space, covarray.sizex);
	for(ordinal i = 0; i < covarray.sizey; i++)
		for(ordinal j = 0; j < covarray.sizex; j++)
		{
			covarray.data[i][j] /= (Data.sizey-1);
			if (j == i) sum += covarray.data[i][j];
		}
	LAPACKE_dsyevr(LAPACK_ROW_MAJOR, 'V', 'I', 'U', covarray.sizex, covarray.space, covarray.sizex, vl, vu, il, iu, -1.0, &numvalues, evalues.data, evectors.space, evectors.sizex, suppz.data);
}

void get_pc(MAT2D<ST> & Data, MAT2D<ST> & evectors, MAT2D<ST> & pcarray)
{
	std::fill(pcarray.space, pcarray.space+(pcarray.sizex*pcarray.sizey), 0.);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, pcarray.sizey, pcarray.sizex, Data.sizex,1, Data.space, Data.sizex, evectors.space, evectors.sizex, 0., pcarray.space, pcarray.sizex);

}

void get_correlation(const MAT2D<ST> & Data, const MAT2D<ST> & PC)
{
	ST sumpc = 0;
	MAT1D<ST> nr(Data.sizex), dr1(Data.sizex);
	std::fill(nr.data, nr.data+nr.sizex, 0.);
	std::fill(dr1.data, dr1.data+dr1.sizex, 0.);
	for(ordinal i = 0; i < Data.sizey; i++)
	{
		sumpc += PC.data[i][1]*PC.data[i][1];
		for(ordinal j = 0; j < Data.sizex; j++)
		{
			nr.data[j] += Data.data[i][j]*PC.data[i][1];
			dr1.data[j] += Data.data[i][j]*Data.data[i][j];
		}
	}

	fstream f("Correlation.txt", ios::out);
	for(ordinal i = 0; i < Data.sizex; i++)
		f << nr.data[i]/(sqrt(dr1.data[i]*sumpc)) << endl;
}


ordinal main(ordinal argc, char *argv[])
{

	ordinal nprocs, rank;	

  	MPI_Init(&argc,&argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        Teuchos::oblackholestream nullstream;
        ostream &out = (rank == 0)?cout:nullstream;

	Teuchos::CommandLineProcessor clp;
	clp.setDocString("This program computes the PCA and Support vectors of given data\n");
	ordinal nx = 1553;
	clp.setOption("nx", &nx, "Number of variables");

	ordinal ny = 272;
	clp.setOption("ny", &ny, "Number of samples");

	clp.recogniseAllOptions(true);
	clp.throwExceptions(false);

	Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = clp.parse(argc, argv);
	if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
	if(parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) return 1;
		
	RCP<TableHelp> draw = rcp(new TableHelp);
	out << draw->bigline() << endl;
	ifstream datafile("data.txt");
	if (datafile.is_open())
	{
		RCP<Time> iotime = TimeMonitor::getNewTimer("Total time");
		TimeMonitor Timer(*iotime);
		MAT2D<ST> darray(nx,ny);
		MAT1D<ordinal> carray(ny);
		MAT1D<ST> sumarray(nx); fill(sumarray.data, sumarray.data+sumarray.sizex, 0.);
		string crap;
		ordinal cnt, tot0, tot1; cnt = tot0 = tot1 = 0;
		getline(datafile,crap);
		while (getline(datafile, crap))
		{
			istringstream icrap(crap);
			for (ordinal i = 0; i < darray.sizex; i++) 
			{
				icrap >> darray.data[cnt][i];
				sumarray.data[i] += darray.data[cnt][i];
			} 
			icrap >> carray.data[cnt];
			(carray.data[cnt] == 0)?tot0++:tot1++;
			cnt++;
		}
		datafile.close();
		out << "Finished reading data. There are " << tot0 << " class 0 samples" << " and " << tot1 << " class 1 samples" << endl;
		histogram2(darray, 200, "Initial_histogram.dat");
		MAT2D<ST> avgarray(nx, ny);	
		get_mean(darray, avgarray, sumarray);	
		histogram2(avgarray, 200, "Combined_histogram.dat");
		MAT2D<ST> c0darray(nx,tot0), c1darray(nx, tot1);
		MAT1D<ST> m0darray(nx), m1darray(nx);
		std::fill(m0darray.data, m0darray.data+m0darray.sizex, 0.);
		std::fill(m1darray.data, m1darray.data+m1darray.sizex, 0.);
		ordinal tmp0, tmp1; tmp0 = tmp1 = 0;
		for(ordinal i = 0; i < darray.sizey; i++)
		{
			if (carray.data[i] == 0)
			{
				for(ordinal j = 0; j < darray.sizex; j++)					
				{
					m0darray.data[j] += avgarray.data[i][j];
					c0darray.data[tmp0][j] = avgarray.data[i][j];
				}
				tmp0++;
			}	
			else
			{
				for(ordinal j = 0; j < darray.sizex; j++)
				{
					m1darray.data[j] += avgarray.data[i][j];
                                        c1darray.data[tmp1][j] = avgarray.data[i][j];
				}
				tmp1++;
			}
		}
		for(ordinal i = 0; i < nx; i++)
		{
			m0darray.data[i] /= tot0;
			m1darray.data[i] /= tot1;
		}
		MAT2D<ST> evectors(2,avgarray.sizex), pcarray(2, avgarray.sizey);
		MAT1D<ST> evalues(avgarray.sizex);
		get_cv_and_eig(avgarray, evalues, evectors);
		get_pc(avgarray, evectors, pcarray);
		
		write_data_asc(avgarray, "centered_data.txt");
		write_data1d_asc(carray, "class");		
		
		write_data_asc(c0darray, "centered_class0.txt");
		write_data_asc(c1darray, "centered_class1.txt");

		write_data_asc(evectors, "Eigen_vaectors.txt");


		MAT2D<ST> avgpcarray(2, avgarray.sizey);
		get_mean2d(pcarray, avgpcarray);

		get_correlation(avgarray, avgpcarray);

		write_data_asc(avgpcarray, "centered_PC");

		ofstream odatafile2("PC.txt");
		odatafile2 << "PC2" << " " << "PC1" << " " << "Class" << endl; 
              	for (ordinal i = 0; i < pcarray.sizey; i++)
		{
			for(ordinal j = 0; j < pcarray.sizex; j++)
                		odatafile2 << pcarray.data[i][j] << " ";
			odatafile2 << carray.data[i];
			odatafile2 << endl;
		}
                odatafile2.close();

		ofstream odatafile3("Variable_means.txt");
                for(ordinal i = 0; i < darray.sizex; i++)
                	odatafile3 << m0darray.data[i] << "  " << m1darray.data[i] << endl;
                odatafile3.close();


		std::vector<vec_space> vars;
		std::vector<ST> clabel; 
		vec_space tdata;
		for(ordinal i = 0; i < avgarray.sizey; i++)
		{
			for(ordinal j = 0; j < 310; j+=1)
				tdata(j) = avgarray.data[i][j*5];
			vars.push_back(tdata);
			(carray.data[i] == 0)?clabel.push_back(1):clabel.push_back(-1);
		}

		const ST nmax = maximum_nu(clabel);
		svm_nu_trainer<ktype> trainer;
		for(ST i = 0.00001; i < 1.; i *=5)
		{
			for(ST j = 0.00001; j < nmax; j *= 5)
			{
				trainer.set_kernel(ktype(i));
				trainer.set_nu(j);
				cross_validate_trainer(trainer, vars, clabel, 3);
			}

		}

		trainer.set_kernel(ktype(0.00125));
		trainer.set_nu(0.15625);
		out << "Cross validation accuracy: "  << cross_validate_trainer(trainer, vars, clabel, 3);
		func learner;
		learner = trainer.train(vars, clabel);
		 ordinal cnt0, cnt1; cnt0 = cnt1 = 0;
                for(ordinal i = 0; i < avgarray.sizey; i++)
                {
                        for(ordinal j = 0; j < 310; j+=1)
                                tdata(j) = avgarray.data[i][j*5+3];
			if ((carray.data[i] == 0) && (carray.data[i] == ((learner(tdata)>0)?0:1))) cnt0++;
			else if ((carray.data[i] == 1) && (carray.data[i] == ((learner(tdata)>0)?0:1))) cnt1++;
		}

		out << (100.*cnt0)/tot0 << " % class 0 samples identified correctly using a basis of 310 vectors" << endl;
		out << (100.*cnt1)/tot1 << " % class 1 samples identified correctly using a basis of 310 vectors" << endl;

		
/*		

		krr_trainer<ktype> trainer;
		trainer.use_classification_loss_for_loo_cv();
		for (ST gamma = 0.000001; gamma <= 1; gamma *= 5)
		{
			trainer.set_kernel(ktype(gamma));
			std::vector<ST> loo_values; 
			trainer.train(vars, clabel, loo_values);
			const ST classification_accuracy = mean_sign_agreement(clabel, loo_values);
			out << "gamma: " << gamma << "     LOO accuracy: " << classification_accuracy << endl;
		}
*/
		get_variance(c0darray, "C0Variance.txt");
		get_variance(c1darray, "C1Variance.txt");
		get_variance(avgarray, "Total_Variance.txt");
		get_coeffv(c0darray, "Coeff_Variance_0.txt");
		get_coeffv(c1darray, "Coeff_Variance_1.txt");
		write_data(darray, "full_image.raw");	
		write_data(avgarray, "Full_image_avg.raw");
		write_data(c0darray, "class0.raw");
		write_data(c1darray, "class1.raw");
	}
	else
		out << "Couldn't open the file" << endl;

	TimeMonitor::summarize(out,true);
	MPI_Finalize();
}
