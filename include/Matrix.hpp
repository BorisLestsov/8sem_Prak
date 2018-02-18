#ifndef PRAK4_MATRIX_HPP
#define PRAK4_MATRIX_HPP


//-----implementation-------

template <typename T>
Matrix<T>::Matrix():
        arr(NULL),
        rows(0),
        cols(0)
{}

template <typename T>
Matrix<T>::Matrix(istream& in, InType in_t) {
    if (!in.good())
        throw string("Could not open input file");

    if (in_t == Normal) {
        in >> rows >> cols;

        arr = new T[rows*cols];

        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j)
                in >> arr[i*cols + j];
    } else if (in_t == Binary) {
        in.read((char*) &rows, sizeof(rows));
        in.read((char*) &cols, sizeof(cols));

        arr = new T[rows*cols];

        in.read((char*) arr, sizeof(arr[0]) * cols * rows);
    } else throw string("Unknown input type for matrix");
}

template <typename T>
Matrix<T>::Matrix(T* buf, size_t rows, size_t cols):
        rows(rows),
        cols(cols)
{
    arr = new T[cols*rows];
    memcpy(arr, buf, sizeof(T) * rows*cols);
}

template<class T>
Matrix<T>::Matrix(size_t rows, size_t cols, bool is_rand): rows(rows), cols(cols) {
    size_t i, j;

    arr = new T[rows*cols];

    if (is_rand) {
        srand(time(0));
        for (i = 0; i < rows; i++)
            for (j = 0; j < cols; ++j)
                arr[i*cols + j] = rand() / (T)RAND_MAX;
    }
}

template<class T>
Matrix<T>::Matrix(size_t rows, size_t cols, T elem): rows(rows), cols(cols) {
    size_t i, j;

    arr = new T[cols*rows];

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            arr[i*cols+j] = elem;
}

template<class T>
Matrix<T>::~Matrix() {
    delete []arr;
    arr = NULL;
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &matr): rows(matr.rows), cols(matr.cols) {
    arr = new T[matr.rows*matr.cols];

    memcpy(arr, matr.arr, sizeof(matr.arr[0])*matr.rows*matr.cols);
}

/*
template<class T>
Matrix<T>::Matrix(Matrix&& m):
        arr(m.arr),
        rows(m.rows),
        cols(m.cols)
{
    m.rows = 0;
    m.cols = 0;
    m.arr = NULL;
}
*/

template<class T>
inline size_t Matrix<T>::n_rows() const{
    return rows;
}

template<class T>
inline size_t Matrix<T>::n_cols() const{
    return cols;
}

template<class T>
inline size_t Matrix<T>::size() const{
    return rows*cols;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix &matr) throw(string) {
    try {
        if (rows != matr.rows || cols != matr.cols) throw string("Matrices of different sizes in \"+\"");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<T> tmp = *this;
    int i, j;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            tmp.arr[i*cols+j] += matr.arr[i*cols+j];
    return tmp;
}


template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &matr) throw(string) {
    try {
        if (rows != matr.rows || cols != matr.cols) throw string("Matrices of different sizes in \"-\"");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<T> tmp = *this;
    int i, j;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            tmp.arr[i*cols + j] -= matr.arr[i*cols + j];
    return tmp;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix &matr) throw(string) {
    try {
        if (cols != matr.rows) throw string("Matrices of wrong sizes in \"-\"");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<T> tmp(rows, matr.cols, 0.0f);
    size_t i, j, k;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            for (k = 0; k < cols; k++) {
                tmp.arr[i * cols + j] += arr[i * cols + k] * matr.arr[k * cols + j];
            }
    return tmp;
}

template<class T>
T &Matrix<T>::operator()(size_t i, size_t j) const throw(string) {
    try {
        if (i >= rows || j >= cols || i < 0 || j < 0) throw string("Out of bounds");
    }
    catch (const string& s) {
        cerr << s << endl;
        exit(-1);
    }
    return arr[i*cols + j];
}

template<class Y>
const istream& operator>>(const istream &strm, Matrix<Y> &matr) {
    int i, j;

    for (i = 0; i < matr.rows; i++)
        for (j = 0; j < matr.cols; j++)
            cin >> matr.arr[i*matr.cols + j];
    return strm;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix &m) {
    size_t i, j;

    delete []arr;

    rows = m.n_rows();
    cols = m.n_cols();
    arr = new T[m.n_rows() * m.n_cols()];

    memcpy(arr, m.data(), m.n_rows()*m.n_cols()*sizeof(T));
    return *this;
}
/*
template<class T>
Matrix<T>& Matrix<T>::operator=(Matrix&& m) {
    delete []arr;

    arr = m.arr;
    rows = m.rows;
    cols = m.cols;

    m.rows = 0;
    m.cols = 0;
    m.arr = NULL;
    return *this;
}
*/
template<typename T>
T* Matrix<T>::data() const {
    return arr;
}

template<typename T>
void Matrix<T>::write(ostream& o) const {
    o.write((char*) &rows, sizeof(rows));
    o.write((char*) &cols, sizeof(cols));

    o.write((char*) arr, sizeof(arr[0]) * cols * rows);
}

template<typename T>
void Matrix<T>::print(OutType o_type,
                      uint precision,
                      ostream &stream) const throw(string) {
    uint width;

    switch (o_type) {
        case O_SCI:
            stream << scientific << setprecision(precision);
            width = precision + 12;
            break;
        case O_FXD:
            stream << fixed << setprecision(precision);
            width = precision + 8;
            break;
        default:
            throw string("Wrong output mode");
    }
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            stream << setw(width) << arr[i * cols + j];
        }
        stream << endl;
    }
    stream << endl;
}

template <typename T>
double Matrix<T>::norm() {
    double sum = 0.0;

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++) {
            sum += std::abs(arr[i*cols + j]);
        }
    return sum;
}



// COMPLEX
// NOT SAFE


template <typename T>
Matrix<std::complex<T> >::Matrix():
        arr(NULL),
        rows(0),
        cols(0)
{}

template <typename T>
Matrix<std::complex<T> >::Matrix(istream& in, InType in_t) {
    if (!in.good())
        throw string("Could not open input file");

    if (in_t == Normal) {
        in >> rows >> cols;

        arr = new std::complex<T>[rows*cols];

        for (size_t i = 0; i < rows; ++i)
            for (size_t j = 0; j < cols; ++j)
                in >> arr[i*cols + j];
    } else if (in_t == Binary) {
        in.read((char*) &rows, sizeof(rows));
        in.read((char*) &cols, sizeof(cols));

        arr = new std::complex<T>[rows*cols];

        in.read((char*) arr, sizeof(arr[0]) * cols * rows);
    } else throw string("Unknown input type for matrix");
}

template <typename T>
Matrix<std::complex<T> >::Matrix(std::complex<T>* buf, size_t rows, size_t cols):
        rows(rows),
        cols(cols)
{
    arr = new std::complex<T>[cols*rows];
    memcpy(arr, buf, sizeof(std::complex<T>) * rows*cols);
}

template<class T>
Matrix< std::complex<T> >::Matrix(size_t rows, size_t cols, bool is_rand): rows(rows), cols(cols) {
    size_t i, j;

    arr = new std::complex<T>[rows*cols];

    if (is_rand) {
        srand(time(0));
        for (i = 0; i < rows; i++)
            for (j = 0; j < cols; ++j) {
                arr[i * cols + j].real(rand() / (T) RAND_MAX);
                arr[i * cols + j].imag(rand() / (T) RAND_MAX);
            }
    }
}

template<class T>
Matrix<std::complex<T> >::Matrix(size_t rows, size_t cols, std::complex<T> elem): rows(rows), cols(cols) {
    size_t i, j;

    arr = new std::complex<T>[cols*rows];

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            arr[i*cols+j] = elem;
}

template<class T>
Matrix<std::complex<T> >::~Matrix() {
    delete []arr;
    arr = NULL;
}

template<class T>
Matrix<std::complex<T> >::Matrix(const Matrix<std::complex<T> > &matr): rows(matr.rows), cols(matr.cols) {
    arr = new std::complex<T>[matr.rows*matr.cols];

    memcpy(arr, matr.arr, sizeof(matr.arr[0])*matr.rows*matr.cols);
}

template<class T>
inline size_t Matrix<std::complex<T> >::n_rows() const{
    return rows;
}

template<class T>
inline size_t Matrix<std::complex<T> >::n_cols() const{
    return cols;
}

template<class T>
inline size_t Matrix<std::complex<T> >::size() const{
    return rows*cols;
}

template<class T>
Matrix<std::complex<T> > Matrix<std::complex<T> >::operator+(const Matrix &matr) throw(string) {
    try {
        if (rows != matr.rows || cols != matr.cols) throw string("Matrices of different sizes in \"+\"");
    }
    catch (const string &s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<T> tmp = *this;
    int i, j;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            tmp.arr[i*cols+j] += matr.arr[i*cols+j];
    return tmp;
}


template<class T>
Matrix<std::complex<T> > Matrix<std::complex<T> >::operator-(const Matrix &matr) throw(string) {
    try {
        if (rows != matr.rows || cols != matr.cols) throw string("Matrices of different sizes in \"-\"");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<T> tmp = *this;
    int i, j;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            tmp.arr[i*cols + j] -= matr.arr[i*cols + j];
    return tmp;
}

template<class T>
Matrix<std::complex<T> > Matrix<std::complex<T> >::operator*(const Matrix &matr) throw(string) {
    try {
        if (cols != matr.rows) throw string("Matrices of wrong sizes in \"-\"");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }

    Matrix<std::complex<T> > tmp(rows, matr.cols, 0.0);
    size_t i, j, k;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            for (k = 0; k < cols; k++) {
                tmp(i, j) += this->operator()(i, k) * matr(k, j);
            }
    return tmp;
}

template<class T>
std::complex<T> &Matrix<std::complex<T> >::operator()(size_t i, size_t j) const throw(string) {
    try {
        if (i >= rows || j >= cols || i < 0 || j < 0) throw string("Out of bounds");
    }
    catch (const string *s) {
        cerr << s << endl;
        exit(-1);
    }
    return arr[i*cols + j];
}

template<class Y>
const istream& operator>>(const istream &strm, Matrix<std::complex<Y> > &matr) {
    int i, j;

    for (i = 0; i < matr.rows; i++)
        for (j = 0; j < matr.cols; j++)
            cin >> matr.arr[i*matr.cols + j];
    return strm;
}

template<class T>
Matrix<std::complex<T> > &Matrix<std::complex<T> >::operator=(const Matrix<std::complex<T> > &m) {
    size_t i, j;

    delete []arr;

    rows = m.n_rows();
    cols = m.n_cols();
    arr = new std::complex<T>[m.n_rows() * m.n_cols()];

    memcpy(arr, m.data(), m.n_rows()*m.n_cols()*sizeof(complex<T>));
    return *this;
}

template<typename T>
std::complex<T>* Matrix<std::complex<T> >::data() const {
    return arr;
}

template<typename T>
void Matrix<std::complex<T> >::write(ostream& o) const {
    o.write((char*) &rows, sizeof(rows));
    o.write((char*) &cols, sizeof(cols));

    o.write((char*) arr, sizeof(arr[0]) * cols * rows);
}

template<typename T>
void Matrix<std::complex<T> >::print(OutType o_type,
                      uint precision,
                      ostream &stream) const throw(string) {
    uint width;

    switch (o_type) {
        case O_SCI:
            stream << scientific << setprecision(precision);
            width = precision + 12;
            break;
        case O_FXD:
            stream << fixed << setprecision(precision);
            width = precision + 8;
            break;
        default:
            throw string("Wrong output mode");
    }
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            stream << setw(width) << arr[i * cols + j];
        }
        stream << endl;
    }
    stream << endl;
}

template <typename T>
double Matrix<std::complex<T> >::norm() {
    double sum = 0.0;

    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++) {
            sum += std::abs(arr[i*cols + j]);
        }
    return sum;
}



#endif //PRAK4_MATRIX_HPP
