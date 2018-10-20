# How to contribute

It is recommended to discuss the change you wish to make via issue or email with the authors of Linopt before contributing. Then follow the standard procedure to create a pull request:

1. Fork the repository on GitHub.
2. Create a topic branch based on the `master` branch.
3. Perform modifications taking into account guidelines provided in the sections below.
4. Commit your modifications.
5. Create a pull request.

## Requirements

* Use C++11 standard. C++14 or C++17 is NOT allowed without extreme necessity.
* The code should be cross platform. Usually this means that it uses only STL and cross-platform libraries and does not include OS specific headers.
* Cross platform implies that code should be cross compiler. In case when you need to include compiler-specific commands, e.g., bultins for optimization, use conditional compilation to select proper commands supported by a target compiler. An example of this approach can be found in [lib/misc.h](lib/misc.h) for a `ctz()` function.
* Ensure that indentation is done via tabulation that is 4 spaces long.

## Documentation

Please provide documentation to new functions, classes, class members, typedefs, global constants, etc. as a commentary block before the entity you comment in a [Doxygen](http://www.doxygen.org/) format. The commentary block should not exceed 80-symbols width.

The following prototype of a commentary block is recommended:
```
/** @ingroup group
 * @brief Brief description.
 *
 * @param[in] p1 -- input parameter.
 * @param[out] p2 -- output parameter.
 * @return Calculated value.
 *
 * Detailed description.
 *
 * @throw
 * Throw `exception_class` in case of error.
 *
 * @todo
 * Some stuff.
 *
 * @see bar()
 *
 * @see
 * A. B. Author "Article title." Journal Title __1__, pp. 1-10 (2018),
 * http://doi.org/xxxx
 */
T foo(const T1 &p1, T2 &p2)
{
    // Function body
}
```
Include only those Doxygen commands that are applicable. If the description is short and the entity you document is self-explanatory, only brief description may be given. You can type it as either:
```
/**
 * @brief Brief description.
 */
```
or using one-line notation:
```
/// Brief description.
```