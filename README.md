# UVD-PSU

Reference implementation of the **updatable and verifiable delegated private set union (UVD-PSU)** protocol from the paper
*Verifiable delegation of private set unions over dynamic sets*
by Myungsun Kim (submitted to IEEE TIFS).

The protocol enables two clients (Alice and Bob) to outsource their private sets to a single server which obliviously computes the union, while
* keeping each client's set private from the server,
* allowing a designated receiver to verify the correctness of the returned union, and
* supporting per-element insertions and deletions on the outsourced sets at $O(k)$ symmetric-key cost.

## Building

### Dependencies

* **OpenFHE** >= 1.5.0 - BFV homomorphic encryption (https://github.com/openfheorg/openfhe-development)
* **NTL** - large-integer arithmetic (PRF index reduction)
* **GMP** - multi-precision arithmetic (transitive dependency of NTL)
* **Crypto++** >= 8.9 - HMAC-SHA256 PRF
* **OpenMP** - multi-threaded encoding/decoding (parallelism is required for the reported numbers)

On macOS:
```bash
brew install ntl gmp cryptopp libomp
```

OpenFHE must be built from source.

### Default (64-bit) build

```bash
mkdir build && cd build
cmake ..
make -j
```

### 128-bit `NATIVE_SIZE` build

For experiments with the 128-bit native-integer OpenFHE build, set `VDPSU_NATIVE128=ON`. The CMake variable `OpenFHE_DIR` is then set automatically to `/Users/$USER/opt/openfhe-128/lib/OpenFHE`; adjust `CMakeLists.txt` if your install path differs.

```bash
mkdir build-128 && cd build-128
cmake .. -DVDPSU_NATIVE128=ON
make -j
```

## Running the benchmark

```bash
./bench -n 1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576 \
        -t 3 -i 0.5 -o results.csv
```

Options:
* `-n` comma-separated list of set sizes to benchmark (defaults to `1000,2000,4000`).
* `-t` number of trials per set size (default 3).
* `-i` intersection ratio `|X_A cap X_B| / n` (default 0.5).
* `-o` CSV output path (default stdout).
* `-v` verbose per-trial output.

The CSV columns are
```
n,encodeA_ms,encodeB_ms,delegateA_ms,delegateB_ms,compute_ms,decode_ms,update_us,total_ms,comm_kb
```
where `update_us` is the mean per-update wall-clock time (Algorithm 5).

## Repository layout

| Path | Contents |
| --- | --- |
| `include/` | Public headers: `params.h`, `client.h`, `cloud.h`, `fhe_utils.h`, `mht.h`, `phf.h`, `prf.h`, `bloom_filter.h`, `timer.h` |
| `src/` | Library implementation (mirrors the algorithms of Section IV of the paper) |
| `bench/main.cc` | Benchmark driver |
| `test/` | Unit tests and parameter probes |
| `CMakeLists.txt` | Build configuration |

## Implementation notes

The prototype uses a **single flat MHT** of size `m = 2n` in place of the paper's `eta` per-bin MHTs. The two are operationally equivalent for one execution (same total slots, same number of ciphertexts, same operation counts); the per-bin structure of Section IV is needed only for the analytical update-lifespan argument and is recovered at runtime by the fact that the flat MHT is naturally packed across `ceil(2n/N)` BFV ciphertexts.

The Bloom filter implementation uses standard insert-only bit arrays. The paper's *Update* algorithm (Algorithm 5) requires deletion, which a Dynamic Bloom Filter (Wang et al., WHY+22) would provide; the present prototype omits BF deletes for simplicity since they do not affect the timing of the dominant work.

## License

Provided for academic use accompanying the paper.

## Contact

Myungsun Kim &mdash; `msunkim@gachon.ac.kr`
