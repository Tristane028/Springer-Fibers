<h1>Springer Fibers — Algebraic Geometry &amp; Combinatorics</h1>

<p><strong>Authors:</strong> Tristan Endo, Alexander Ryan<br>
<strong>Advisor:</strong> Raymond O. Chou<br>
<strong>Institution:</strong> University of California, San Diego<br>
<strong>Year:</strong> 2026 — Present</p>

<blockquote>
  ⚠️ This repository is actively being updated as new research is conducted. Results, implementations, and documentation are subject to change.
</blockquote>

<h2>Overview</h2>

<p>This project investigates the geometry and combinatorics of Springer fibers. We examine the structure of flag varieties under nilpotent operators defined by Jordan blocks, compute valid Young tableau fillings, and analyze the corresponding Schubert cells. The goal of this research is to characterize the irreducible components of Springer fibers through Plücker coordinates, Springer span checks, closure relations, and explicit defining equations.</p>

<h2>What the Scripts Do</h2>

<p>Each script takes a number string as input (e.g. <code>"221"</code>) where each digit defines the size of a Jordan block. The codebase is split into four scripts, each building on the previous one.</p>

<h3>1. <code>Combo Research.py</code> — Schubert cells &amp; Plücker coordinates</h3>
<ul>
  <li>Constructs the Young tableau and Jordan matrix from the input</li>
  <li>Generates all valid fillings of the tableau respecting the order-preserving condition within each Jordan block</li>
  <li>Converts each filling into a Schubert cell matrix with symbolic entries</li>
  <li>Runs Springer span checks to determine free variables and constrained relations</li>
  <li>Computes Plücker coordinates for each cell using minors of the matrix</li>
  <li>Outputs a full summary table of all cells with their tableau, matrix, free variables, relations, and Plücker coordinates</li>
</ul>

<h3>2. <code>Closure Representation.py</code> — Closure relations between cells</h3>
<ul>
  <li>For each pair of cells, determines whether one is contained in the closure of the other via parameter substitutions and limit operations</li>
  <li>Searches over polynomial substitutions (<code>a_i → c · v^k</code> for k up to 3), target parameter values, and a limit variable sent to ∞</li>
  <li>Uses a ratio-limit check to verify projective convergence of Plücker coordinates</li>
  <li>Produces three outputs:
    <ul>
      <li><strong>Closure table</strong> — grid showing which cells lie in the closure of which</li>
      <li><strong>Closure details</strong> — for each closure relation, the exact substitution, limit, and scaled Plücker vectors that realize it</li>
      <li><strong>Closure summary</strong> — for each permutation, the list of cells it contains in its closure</li>
    </ul>
  </li>
</ul>

<h3>3. <code>Locus Construction.py</code> — Parameter spaces &amp; boundary stratification</h3>
<ul>
  <li>For each cell, presents its parameter space 𝔸<sup>d</sup> explicitly with its free coordinates</li>
  <li>Displays the Plücker map φ<sub>w</sub> : 𝔸<sup>d</sup> → 𝔾(•, n) sending free parameters to the Plücker coordinates</li>
  <li>Describes the boundary stratification: for each cell in the closure, gives the locus inside 𝔸<sup>d</sup> (substitutions + limit direction) and the target parameters at the limit point</li>
</ul>

<h3>4. <code>Plucker Relations.py</code> — Plücker relations &amp; cell ideals</h3>
<ul>
  <li>Generates the <strong>universal Plücker relations</strong> for the relevant Grassmannian 𝔾(k, n) (one quadratic relation per (k+2)-subset)</li>
  <li>Generates <strong>incidence relations</strong> between consecutive levels (e.g. 1-Plückers to 2-Plückers) coming from the flag condition V<sub>r</sub> ⊂ V<sub>r+1</sub></li>
  <li>Verifies that each cell's computed Plücker coordinates satisfy every universal relation</li>
  <li>Computes the <strong>cell ideal</strong> for each cell via Gröbner-basis elimination: the polynomial relations among Plücker coordinates that cut out the cell's image in the Grassmannian. Most relations are linear (<code>p_I = 0</code> or <code>p_I = ±1</code>), but higher-dimensional cells produce genuine polynomial relations (e.g. <code>p2_12 = (p3_123)²</code> for the top cell of (12)(34))</li>
</ul>

<h2>Running the Scripts</h2>

<pre><code>python "Combo Research.py"
python "Closure Representation.py"
python "Locus Construction.py"
python "Plucker Relations.py"
</code></pre>

<p>Each script will prompt for a number string:</p>

<pre><code>Enter number string: 221
</code></pre>

<p>The scripts can be run independently; each carries out the cell-generation and Plücker-coordinate computation on its own before adding its specialized analysis.</p>

<h2>Dependencies</h2>
<ul>
  <li>Python 3.9</li>
  <li><code>numpy</code></li>
  <li><code>sympy</code></li>
</ul>

<p>Install with:</p>

<pre><code>pip install numpy sympy
</code></pre>

<h2>Collaborators</h2>
<ul>
  <li><a href="https://github.com/raymondchou420">Raymond Chou</a></li>
  <li><a href="https://github.com/Tristane028">Tristan Endo</a></li>
  <li><a href="https://github.com/Wave449">Alexander Ryan</a></li>
</ul>

<h2>References</h2>
<p>Talia Goldwasser, Meera Nadeem, Garcia Sun, and Julianna Tymoczko. <em>Cell Closures for Two-Row Springer Fibers via Noncrossing Matchings.</em> 📄 <a href="https://docs.google.com/viewer?url=https://raw.githubusercontent.com/Tristane028/Springer-Fibers/main/Cell_Closures_for_Two-Row_Springer_Fibers.pdf">View Paper</a></p>
