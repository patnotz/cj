1. Create a build directory
    <pre><code>
    $ mkdir build
    $ cd build
    </code></pre>
1. Copy and tailor your CMake script
    <pre><code>
    $ cp ../doc/do-cmake
    $ emacs do-cmake # set your trilinos installation location
    </code></pre>
1. Build
    <pre><code>
    $ ./do-cmake
    $ make
    </code></pre>
1. Run it
    <pre><code>
    src/cj
    </code></pre>
