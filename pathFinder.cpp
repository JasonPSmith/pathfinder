//INPUTS:
// vertices: address of numpy file of indexes of vertices to be considered (will consider all 2 element subsets) dtype uint32
// threads: number of parallel threads
// size: number of vertices in graph
// pathdist: Address of numpy file containing path distance between the elements to consider should be of size |vertices|^2 dtype uint32
// indices: Address of numpy file of indices of csr format of graph dtype uint32
// indptr: Address of numpy file of indptr of csr format of graph dtype uint32
// out: Output address will save a file for each thread

// LIBRARIES
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <vector>
#include <thread>
#include <functional>
#include <deque>
#include <sstream>
#include <cmath>
#include <chrono>

// Libraries needed for cnpy
#include <stdexcept>
#include <cstdio>
#include <typeinfo>
#include <cassert>
#include <zlib.h>
#include <map>
#include <memory>
#include <stdint.h>
#include <numeric>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <stdexcept>
#include <regex>


//##############################################################################
// TYPE DEFINTIONS
typedef uint64_t vertex_t;//Requires numpy arrays to be given as uint64 dtype


//##############################################################################
// Class for thread safe printing
class LockIO {
    static pthread_mutex_t *mutex;
    public:
        LockIO() { pthread_mutex_lock( mutex ); }
        ~LockIO() { pthread_mutex_unlock( mutex ); }
};

static pthread_mutex_t* getMutex() {
    pthread_mutex_t *mutex = new pthread_mutex_t;
    pthread_mutex_init( mutex, NULL );
    return mutex;
}
pthread_mutex_t* LockIO::mutex = getMutex();


//********************************************************************//

struct NpyArray {
    NpyArray(const std::vector<size_t>& _shape, size_t _word_size, bool _fortran_order) :
        shape(_shape), word_size(_word_size), fortran_order(_fortran_order)
    {
        num_vals = 1;
        for(size_t i = 0;i < shape.size();i++) num_vals *= shape[i];
        data_holder = std::shared_ptr<std::vector<char>>(
            new std::vector<char>(num_vals * word_size));
    }

    NpyArray() : shape(0), word_size(0), fortran_order(0), num_vals(0) { }

    template<typename T>
    T* data() {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    const T* data() const {
        return reinterpret_cast<T*>(&(*data_holder)[0]);
    }

    template<typename T>
    std::vector<T> as_vec() const {
        const T* p = data<T>();
        return std::vector<T>(p, p+num_vals);
    }

    size_t num_bytes() const {
        return data_holder->size();
    }

    std::shared_ptr<std::vector<char>> data_holder;
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    size_t num_vals;
};

void parse_npy_header(FILE* fp,size_t& word_size, std::vector<size_t>& shape, bool& fortran_order);
void parse_npy_header(unsigned char* buffer,size_t& word_size, std::vector<size_t>& shape, bool& fortran_order);
NpyArray npy_load(std::string fname);

void parse_npy_header(unsigned char* buffer,size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {
    //std::string magic_string(buffer,6);
    uint8_t major_version = *reinterpret_cast<uint8_t*>(buffer+6);
    uint8_t minor_version = *reinterpret_cast<uint8_t*>(buffer+7);
    uint16_t header_len = *reinterpret_cast<uint16_t*>(buffer+8);
    std::string header(reinterpret_cast<char*>(buffer+9),header_len);

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order")+16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stol(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable.
    //not sure when this applies except for byte array
    loc1 = header.find("descr")+9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

void parse_npy_header(FILE* fp, size_t& word_size, std::vector<size_t>& shape, bool& fortran_order) {
    char buffer[256];
    size_t res = fread(buffer,sizeof(char),11,fp);
    if(res != 11)
        throw std::runtime_error("parse_npy_header: failed fread");
    std::string header = fgets(buffer,256,fp);
    assert(header[header.size()-1] == '\n');

    size_t loc1, loc2;

    //fortran order
    loc1 = header.find("fortran_order");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'fortran_order'");
    loc1 += 16;
    fortran_order = (header.substr(loc1,4) == "True" ? true : false);

    //shape
    loc1 = header.find("(");
    loc2 = header.find(")");
    if (loc1 == std::string::npos || loc2 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: '(' or ')'");

    std::regex num_regex("[0-9][0-9]*");
    std::smatch sm;
    shape.clear();

    std::string str_shape = header.substr(loc1+1,loc2-loc1-1);
    while(std::regex_search(str_shape, sm, num_regex)) {
        shape.push_back(std::stol(sm[0].str()));
        str_shape = sm.suffix().str();
    }

    //endian, word size, data type
    //byte order code | stands for not applicable.
    //not sure when this applies except for byte array
    loc1 = header.find("descr");
    if (loc1 == std::string::npos)
        throw std::runtime_error("parse_npy_header: failed to find header keyword: 'descr'");
    loc1 += 9;
    bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
    assert(littleEndian);

    //char type = header[loc1+1];
    //assert(type == map_type(T));

    std::string str_ws = header.substr(loc1+2);
    loc2 = str_ws.find("'");
    word_size = atoi(str_ws.substr(0,loc2).c_str());
}

NpyArray load_the_npy_file(FILE* fp) {
    std::vector<size_t> shape;
    size_t word_size;
    bool fortran_order;
    parse_npy_header(fp,word_size,shape,fortran_order);

    NpyArray arr(shape, word_size, fortran_order);
    size_t nread = fread(arr.data<char>(),1,arr.num_bytes(),fp);
    if(nread != arr.num_bytes())
        throw std::runtime_error("load_the_npy_file: failed fread");
    return arr;
}


NpyArray npy_load(std::string fname) {

    FILE* fp = fopen(fname.c_str(), "rb");

    if(!fp) throw std::runtime_error("npy_load: Unable to open file "+fname);

    NpyArray arr = load_the_npy_file(fp);

    fclose(fp);
    return arr;
}

//******************************************************************//
//Input Functions

typedef std::unordered_map<std::string, const char*> named_arguments_t;

//Parse the inputted argument by creating pairs (arg_name, value)
const named_arguments_t parse_arguments(int argc, char** argv) {
    named_arguments_t named_arguments;

    for (size_t i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        arg = arg.substr(2); // remove --
        named_arguments.insert(std::make_pair(arg, argv[++i]));  // We have extra data, so parse it
    }

    return named_arguments;
}

//A structure that contains all the inputted data, and stores the results
struct parameters_t {
    unsigned short parallel_threads = 8;
    vertex_t number_of_vertices = 0;
    vertex_t vertices_completed = 0;
    std::vector<vertex_t> vertices;
    std::vector<vertex_t> maxpath;
    std::vector<vertex_t> indices;
    std::vector<vertex_t> indptr;
    std::vector<std::vector<vertex_t>> path_dist;
    std::vector<std::ofstream> outstreams;


    //Constructor
    parameters_t(int argc, char** argv){
        named_arguments_t named_arguments = parse_arguments(argc, argv);

        //iterators to find inputted arguments
        auto it_vertices = named_arguments.find("vertices");
        auto it_threads = named_arguments.find("threads");
        auto it_size = named_arguments.find("size");
        auto it_pathdist = named_arguments.find("pathdist");
        auto it_indices = named_arguments.find("indices");
        auto it_indptr = named_arguments.find("indptr");
        auto it_out = named_arguments.find("out");
        named_arguments_t::const_iterator it;

        number_of_vertices = atol(it_size->second);
        parallel_threads = atoi(it_threads->second);
        if (it_vertices == named_arguments.end() || it_threads == named_arguments.end() || it_size == named_arguments.end() || it_pathdist == named_arguments.end() || it_out == named_arguments.end() || it_indices == named_arguments.end() || it_indptr == named_arguments.end()) {
            std::cerr << "ERROR: an input is missing" << std::endl; exit(-1);
        }

        NpyArray indices_file = npy_load(it_indices->second);
        NpyArray indptr_file = npy_load(it_indptr->second);
        NpyArray pathdist_file = npy_load(it_pathdist->second);
        NpyArray vertices_file = npy_load(it_vertices->second);
        indices = indices_file.as_vec<vertex_t>();
        indptr = indptr_file.as_vec<vertex_t>();
        vertices = vertices_file.as_vec<vertex_t>();
        std::vector<vertex_t> pathdist_flat = pathdist_file.as_vec<vertex_t>();
        path_dist.resize(vertices.size());
        maxpath.resize(vertices.size());
        for(int i = 0; i < pathdist_flat.size(); i++){
             path_dist[i/vertices.size()].push_back(pathdist_flat[i]);
        }
        for(int i = 0; i < vertices.size(); i++){
             maxpath[i] = *std::max_element(path_dist[i].begin(),path_dist[i].end());
        }

        for (int i = 0; i < parallel_threads; i++){
            outstreams.push_back(std::ofstream(it_out->second+std::to_string(i)+".paths"));
        }
    } //end constructor
};

//******************************************************************//
//Graph Functions

bool is_connected_by_an_edge(vertex_t from, vertex_t to, parameters_t& params) {
    return std::find(params.indices.begin()+params.indptr[from],params.indices.begin()+params.indptr[from+1],to) != params.indices.begin()+params.indptr[from+1];
}

vertex_t get_outgoing_start(vertex_t from, parameters_t& params){
    return params.indptr[from];
}
vertex_t get_outgoing_end(vertex_t from, parameters_t& params){
    return params.indptr[from+1];
}


//*********************************************************************//
//Path Counting Functions


//Checks if array A contains element a
bool array_contains(vertex_t a, std::vector<vertex_t>& A){
    return std::find(std::begin(A), std::end(A), a) != std::end(A);
}
int array_index(vertex_t a, std::vector<vertex_t>& A){
    auto it = find(A.begin(), A.end(), a);
    if (it != A.end()) { return std::distance(A.begin(), it); }
    else { return -1; }
}

//iterator that adds next vertex to possible cycle
void extend_cycle(vertex_t vertex, std::vector<vertex_t>& prefix,
                       int thread_id, parameters_t& params, int loc_start) {

    // Initialise new chords and prefix for this branch
    std::vector<vertex_t> current_prefix(prefix);

    // Add current vertex to prefix
    current_prefix.push_back(vertex);

    if(prefix.size() > 0){

          int loc = array_index(vertex,params.vertices);

          if(loc != -1 && prefix.size() <= params.path_dist[loc_start][loc]){
                for (auto i : current_prefix) params.outstreams[thread_id] << i << " ";
                params.outstreams[thread_id] << std::endl;
          }
    }
    if (current_prefix.size() <= params.maxpath[loc_start]){
        // Recurse through all remaining out_neighs and extend path
        for(auto iter = params.indptr[vertex]; iter != params.indptr[vertex+1]; ++iter){
            if(!array_contains(params.indices[iter],current_prefix)){
                extend_cycle(params.indices[iter], current_prefix, thread_id, params, loc_start);
            }
        }
    }
}

//Assign to each thread the vertices to be considered as source and start the thread computing
void start_count(int thread_id, parameters_t& params) {

    //To this thread assign all vertices v such that v (mod threads) = thread_id
    std::vector<vertex_t> first_vertices;
    for (size_t index = thread_id; index < params.vertices.size(); index += params.parallel_threads){
        first_vertices.push_back(index);
    }

    std::vector<vertex_t> prefix; //vertices in path so far

    //For each element index of first_vertices count the cycles where index is the smallest vertex in the cycle
    for (auto index : first_vertices){
            extend_cycle(params.vertices[index], prefix, thread_id, params, index);
        {
            LockIO lock;
            params.vertices_completed++;
            std::cout << index << " : thread " << thread_id << " : ";
            std::cout << params.vertices_completed << "/" << params.vertices.size() << std::endl;
        }
    }
}


int main(int argc, char** argv) {
    //Parse input
    std::cout << "Reading input...";
    parameters_t params(argc, argv);

    std::cout << "Done" << std::endl << "Counting Paths" << std::endl;

    //Compute the counts
    //Create a vector of threads
    std::vector<std::thread> t(params.parallel_threads - 1);

    //Start each thread counting paths
    for (size_t index = 0; index < params.parallel_threads - 1; ++index){
        t[index] = std::thread(&start_count, index, std::ref(params));
    }
    start_count(params.parallel_threads - 1, params); // Also do work in this thread

    // Wait until all threads stop
    for (size_t i = 0; i < params.parallel_threads - 1; ++i) t[i].join();

    std::cout << "Done" << std::endl;
}
