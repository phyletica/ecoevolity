#ifndef ECOEVOLITY_DEBUG_HPP
#define ECOEVOLITY_DEBUG_HPP

#if defined(IGNORE_ECOEVOLITY_DEBUG) || defined(NDEBUG)
#   define ECOEVOLITY_DEBUG(x)
#else
#   define ECOEVOLITY_DEBUG(x) x
#endif

#endif

