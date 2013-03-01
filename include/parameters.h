#ifndef PARAMETERS_H_
#define PARAMETERS_H_

template<typename T>
struct Parameters
{
};

template<typename T>
void operator<<(Parameters<T>& p, const T& v);

template<typename T>
void operator>>(const Parameters<T>& p, T& v);

template<typename SourceType, typename LensType>
struct SourceLensParameters
{
	Parameters<SourceType> source;
	Parameters<LensType> lens;
};

#endif
