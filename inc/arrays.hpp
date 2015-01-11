using namespace std;
template <class T>
class MAT1D {
		float sizeofdata;
	public:
		T *data;
		MAT1D(long NX);
		void create(long NX);
		float size() const {return sizeofdata;};
		~MAT1D();
		void clean() const;
		long sizex;
};

template <class T>
MAT1D <T>::MAT1D (long NX) {
	if ((data = (T *) mkl_malloc(NX *sizeof(T), 128)) == NULL) 
	{
		cout << "Total memory requested is " << (1.*NX*sizeof(T)/(1024*1024*1024)) << " GB" << endl;
		cerr << "Memory not alloc'd for vector" << endl;
		exit(1);
	}
	sizex = NX;
};

template <class T>
MAT1D <T>::~MAT1D() {
	mkl_free(data);
	data = NULL;
};

template <class T>
void MAT1D <T>::create (long NX) {
        if ((data = (T *) mkl_malloc(NX *sizeof(T), 128)) == NULL) 
	{
                cout << "Total memory requested is " << (1.*NX*sizeof(T)/(1024*1024*1024)) << " GB" << endl;
                cerr << "Memory not alloc'd for vector" << endl;
                exit(1);
        }
      	sizex = NX;
};

template <class T>
void MAT1D <T>::clean() const {
        mkl_free(data);
        data = NULL;
};

template <class T>
class MAT2D {
                float sizeofdata;
        public:
		MAT2D(long NX, long NY);
		~MAT2D(void);
                T **data, *space;
		long sizex, sizey;
                void create(long NX, long NY);
                float size() const {return sizeofdata;};
                void clean() const;           
};

template <class T>
MAT2D<T>::MAT2D (long NX, long NY) {
int y;

        if ((space = (T *)mkl_malloc((NX*NY) * sizeof(T), 128)) == NULL) {
                cout << "Total memory requested is " << (1.*NX*NY*sizeof(T)/(1024*1024*1024)) << " GB" << endl;
                cerr << "Memory not allocd for whole block" << endl;
                exit(1);
        }

        if ((data = (T **)mkl_malloc(NY * sizeof(T *), 128)) == NULL) {
                cerr << "Memory not allocd for y depth wise" << endl;
                exit(1);
        }
        for (y = 0; y < NY; y++) {
                (data)[y] = space + (y * NX);
        }

        sizeofdata = 1.*NX*NY*sizeof(T)/(1024*1024*1024);
	sizex = NX; sizey = NY;
}

template <class T>
MAT2D<T>::~MAT2D() {
        mkl_free(data);
        data = NULL;
        mkl_free(space);
        space = NULL;
}

template <class T>
void MAT2D<T>::create (long NX, long NY) {
int y;

        if ((space = (T *)mkl_malloc((NX*NY) * sizeof(T), 128)) == NULL) {
                cout << "Total memory requested is " << (1.*NX*NY*sizeof(T)/(1024*1024*1024)) << " GB" << endl;
                cerr << "Memory not allocd for whole block" << endl;
                exit(1);
        }

        if ((data = (T **)mkl_malloc(NY * sizeof(T *), 128)) == NULL) {
                cerr << "Memory not allocd for y depth wise" << endl;
                exit(1);
        }
        for (y = 0; y < NY; y++) {
                (data)[y] = space + (y * NX);
        }

        sizeofdata = 1.*NX*NY*sizeof(T)/(1024*1024*1024);
        sizex = NX; sizey = NY;
}

template <class T>
void MAT2D<T>::clean() const {
        mkl_free(data);
        data = NULL;
        mkl_free(space);
        space = NULL;
}
