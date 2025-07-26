export module oml.FakeDouble;
//
//  Make some fake, no-op for double and double containers, so we can use the same template code
//  for complex and double types
//
export {
inline const double& conj(const double& d) { return d;}
inline const double& real(const double& d) { return d;}
inline       double  imag(const double& d) { return 0.0;}
inline       double  Max (const double& d) { return d;}
}