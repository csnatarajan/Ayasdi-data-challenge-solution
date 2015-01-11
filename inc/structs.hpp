struct PComm 
{
MPI_Comm comm;
int rank;
int nprocs;
int coords[3];
long readcount;
bool exist;
long total;
long endsize;
long np;
};
struct velocity
{
float hi;
float lo;
float fact;
float alpha_I;
float g_0;
float alpha_G;
};
