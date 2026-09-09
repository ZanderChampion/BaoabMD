#include <napi.h>
#include <CL/cl.h>
#include <cmath>
#include <vector>

static cl_platform_id platform = NULL;
static cl_device_id device = NULL;
static cl_context context = NULL;
static cl_command_queue queue = NULL;
static cl_program program = NULL;
static cl_kernel kernel = NULL;
static bool initialized = false;

static int current_N = 0;
static cl_mem d_p = NULL, d_f = NULL, d_q = NULL, d_sig = NULL, d_eps = NULL, d_excl = NULL, d_lj = NULL, d_coul = NULL;
static std::vector<float> lj_out;
static std::vector<float> coul_out;

const char* kernelSource = R"CL(
__kernel void compute_forces(
    __global const float* p,
    __global float* f,
    __global const float* q,
    __global const float* sig,
    __global const float* eps,
    __global const unsigned char* excl,
    const float coulombConst,
    const int N,
    __global float* lj_out,
    __global float* coul_out
) {
    int i = get_global_id(0);
    if (i >= N) return;

    int i3 = i * 3;
    float p1x = p[i3]; float p1y = p[i3 + 1]; float p1z = p[i3 + 2];
    float q1 = q[i]; float sig1 = sig[i]; float eps1 = eps[i];

    float fx = 0.0f, fy = 0.0f, fz = 0.0f;
    float local_lj = 0.0f, local_coul = 0.0f;
    
    const float cutoff2 = 100.0f; 

    for (int j = 0; j < N; j++) {
        if (i == j) continue;
        
        unsigned char ex = excl[i * N + j];
        if (ex == 1) continue; 

        int j3 = j * 3;
        float dx = p1x - p[j3];
        float dy = p1y - p[j3 + 1];
        float dz = p1z - p[j3 + 2];
        float r2 = dx * dx + dy * dy + dz * dz + 1e-6f;
        
        if (r2 > cutoff2) continue;

        float r = sqrt(r2);
        float sig_mix = 0.5f * (sig1 + sig[j]);
        float eps_mix = sqrt(eps1 * eps[j]);

        float sr = sig_mix / r;
        float sr2 = sr * sr;
        float sr6 = sr2 * sr2 * sr2;
        float sr12 = sr6 * sr6;
        
        float ljScale = (ex == 2) ? 0.5f : 1.0f;
        float coulScale = (ex == 2) ? 0.833333f : 1.0f;

        float fLJ = ljScale * (24.0f * eps_mix * (2.0f * sr12 - sr6)) / r2;
        local_lj += ljScale * 0.5f * (4.0f * eps_mix * (sr12 - sr6));

        float fCoul = 0.0f;
        if (q1 != 0.0f && q[j] != 0.0f) {
            fCoul = coulScale * (coulombConst * q1 * q[j]) / (r2 * r);
            local_coul += coulScale * 0.5f * ((coulombConst * q1 * q[j]) / r);
        }

        float fNet = fLJ + fCoul;
        fx += dx * fNet;
        fy += dy * fNet;
        fz += dz * fNet;
    }

    f[i3] = fx;
    f[i3 + 1] = fy;
    f[i3 + 2] = fz;
    lj_out[i] = local_lj;
    coul_out[i] = local_coul;
}
)CL";

void CleanBuffers() {
    if (d_p) { clReleaseMemObject(d_p); d_p = NULL; }
    if (d_f) { clReleaseMemObject(d_f); d_f = NULL; }
    if (d_q) { clReleaseMemObject(d_q); d_q = NULL; }
    if (d_sig) { clReleaseMemObject(d_sig); d_sig = NULL; }
    if (d_eps) { clReleaseMemObject(d_eps); d_eps = NULL; }
    if (d_excl) { clReleaseMemObject(d_excl); d_excl = NULL; }
    if (d_lj) { clReleaseMemObject(d_lj); d_lj = NULL; }
    if (d_coul) { clReleaseMemObject(d_coul); d_coul = NULL; }
}

void InitOpenCL() {
    if (initialized) return;
    cl_int err;
    err = clGetPlatformIDs(1, &platform, NULL);
    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    queue = clCreateCommandQueue(context, device, 0, &err);
    program = clCreateProgramWithSource(context, 1, &kernelSource, NULL, &err);
    clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    kernel = clCreateKernel(program, "compute_forces", &err);
    initialized = true;
}

Napi::Value InitStaticData(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    InitOpenCL();

    float* q = info[0].As<Napi::Float32Array>().Data();
    float* sig = info[1].As<Napi::Float32Array>().Data();
    float* eps = info[2].As<Napi::Float32Array>().Data();
    uint8_t* excl = info[3].As<Napi::Uint8Array>().Data();
    size_t N = info[0].As<Napi::Float32Array>().ElementLength();
    int N_int = static_cast<int>(N);

    CleanBuffers();
    current_N = N_int;
    lj_out.resize(N);
    coul_out.resize(N);

    cl_int err;
    d_p = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * N * 3, NULL, &err);
    d_f = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * N * 3, NULL, &err);
    d_q = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N, q, &err);
    d_sig = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N, sig, &err);
    d_eps = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * N, eps, &err);
    d_excl = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(unsigned char) * N * N, excl, &err);
    d_lj = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * N, NULL, &err);
    d_coul = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * N, NULL, &err);

    return env.Null();
}

Napi::Value ComputeForces(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    float* p = info[0].As<Napi::Float32Array>().Data();
    float* f = info[1].As<Napi::Float32Array>().Data();
    float coulombConst = static_cast<float>(info[2].As<Napi::Number>().DoubleValue());

    int N = current_N;
    if (N == 0 || !d_p) return env.Null();

    clEnqueueWriteBuffer(queue, d_p, CL_FALSE, 0, sizeof(float) * N * 3, p, 0, NULL, NULL);

    clSetKernelArg(kernel, 0, sizeof(cl_mem), &d_p);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_f);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_q);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_sig);
    clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_eps);
    clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_excl);
    clSetKernelArg(kernel, 6, sizeof(float), &coulombConst);
    clSetKernelArg(kernel, 7, sizeof(int), &N);
    clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_lj);
    clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_coul);

    size_t globalSize = N;
    clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &globalSize, NULL, 0, NULL, NULL);

    clEnqueueReadBuffer(queue, d_f, CL_FALSE, 0, sizeof(float) * N * 3, f, 0, NULL, NULL);
    clEnqueueReadBuffer(queue, d_lj, CL_FALSE, 0, sizeof(float) * N, lj_out.data(), 0, NULL, NULL);
    clEnqueueReadBuffer(queue, d_coul, CL_TRUE, 0, sizeof(float) * N, coul_out.data(), 0, NULL, NULL);

    double totalLj = 0.0, totalCoulomb = 0.0;
    for (int i = 0; i < N; i++) {
        totalLj += lj_out[i];
        totalCoulomb += coul_out[i];
    }

    Napi::Object result = Napi::Object::New(env);
    result.Set("lj", Napi::Number::New(env, totalLj));
    result.Set("coulomb", Napi::Number::New(env, totalCoulomb));
    return result;
}

Napi::Object Init(Napi::Env env, Napi::Object exports) {
    exports.Set(Napi::String::New(env, "initStaticData"), Napi::Function::New(env, InitStaticData));
    exports.Set(Napi::String::New(env, "computeForces"), Napi::Function::New(env, ComputeForces));
    return exports;
}

NODE_API_MODULE(forces, Init)