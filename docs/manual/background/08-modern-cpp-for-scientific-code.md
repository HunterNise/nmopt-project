# 08 · Modern C++ for scientific code

**Background navigation:** [Index](README.md) \
Previous: [07 · PDE-constrained optimization](07-pde-constrained-optimization.md) \
Next: [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) \
Practical companion: [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md)

## Purpose

Scientific C++ often looks much more complicated than the mathematics it implements. A finite-element program may contain a familiar residual or linear solve, but the surrounding code uses templates, references, smart pointers, virtual interfaces, lambdas, callbacks, and library-specific types. If those language mechanisms are unfamiliar, it is easy to mistake C++ syntax for mathematical complexity.

This chapter develops the subset of modern C++ needed to **read and reason about** interfaces such as those used by `nmopt` and deal.II. The objective is not to turn the reader into a general C++ language expert. Instead, we will learn to answer practical questions while reading scientific code:

- What type of object is being passed here?
- Is it copied, borrowed, moved, or shared?
- Who must keep it alive?
- May this function modify it?
- Is a decision made at compile time through a template or at run time through a virtual interface?
- What does a lambda capture, and can the captured object outlive the lambda?
- Why can several unrelated callable objects all be stored behind one callback type?

The project currently compiles its public C++ interface as C++17. Some references suggested at the end discuss newer C++ standards; the concepts used here are restricted to facilities available in C++17 unless stated otherwise.

This chapter assumes basic programming familiarity: variables, functions, conditionals, loops, and the idea of a user-defined data type. It does not assume previous knowledge of C++ ownership conventions, templates, inheritance, or move semantics.

We deliberately stop before advanced template metaprogramming, allocators, concurrency, low-level memory management, and the full standard library. [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) treats the command line, compilation, linking, Make/Ninja, CMake, and tests separately. [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) then applies the language concepts to a concrete scientific library.

## Before you start

You should be able to read simple procedural code in some language. For example, this should not be mysterious:

```cpp
double squared_norm(const std::vector<double> &x)
{
  double value = 0.0;
  for (const double entry : x)
    value += entry * entry;
  return value;
}
```

You do **not** need to know yet why the parameter is written as `const std::vector<double> &x`; that is one of the first questions we will answer.

The mathematics in this chapter is intentionally light. The focus is on how scientific software expresses mathematical objects and operations.

## What you will be able to do

After this chapter, you should be able to:

- read a C++ function declaration from the outside in;
- distinguish values, references, pointers, and smart pointers;
- interpret `const` on objects, references, member functions, and pointers;
- explain object lifetime and the RAII model at a practical level;
- distinguish ownership from borrowing;
- explain the roles of `std::unique_ptr`, `std::shared_ptr`, and `std::weak_ptr`;
- understand why `std::move` appears when an object is transferred into another object;
- read constructors and member-initializer lists;
- read class templates, function templates, type aliases, and dependent type names;
- distinguish compile-time polymorphism through templates from run-time polymorphism through virtual functions;
- recognize abstract interfaces, pure virtual functions, `override`, and `final`;
- read lambdas and reason about their captures and lifetime requirements;
- understand `std::function` as a type-erased callable wrapper;
- recognize `enum class`, `std::optional`, and `std::variant` in interface code;
- understand the role of exceptions and explicit run-time checks at API boundaries;
- read a callback-based scientific interface end to end;
- recognize the same patterns when they appear in `nmopt` and deal.II code.

## Roadmap

We start with declarations because they are the grammar of an interface: a declaration tells us the name of an object or operation, its type, how values enter, and what comes back. We then separate three concepts that are often conflated by new C++ readers: object identity, lifetime, and ownership.

With those foundations in place, we introduce RAII and smart pointers, followed by move semantics. The second half of the chapter develops the two main abstraction mechanisms that dominate scientific C++: templates, which select types and operations at compile time, and virtual interfaces, which select implementations at run time. Lambdas and `std::function` then show how behavior itself can be passed as data.

We finish with a small callback-based model interface and a short map from these language mechanisms to the current `nmopt` codebase.

## 1. Start by reading the declaration, not the function body

Consider the declaration

```cpp
double energy(const Vector &state);
```

A useful reading order is:

1. the function is named `energy`;
2. it returns a `double`;
3. it receives one argument called `state`;
4. the argument is a reference to a `Vector`;
5. the referred-to vector is `const`, so this interface promises not to modify it through this reference.

Compare it with

```cpp
void scale(Vector &state, double factor);
```

Now the first parameter is a non-const reference. The caller should expect `state` to be modified.

This style of reading scales well. A declaration such as

```cpp
std::shared_ptr<const BlockLayout>
make_layout(const std::vector<std::size_t> &dimensions);
```

may look dense, but it can still be decomposed:

- `make_layout` is a function;
- it receives a read-only reference to a vector of dimensions;
- it returns a `std::shared_ptr`;
- the shared object is a `const BlockLayout`.

The syntax is compact because a declaration carries several independent facts at once.

### 1.1 Types are part of the program's reasoning

In a dynamically typed language one may think primarily about values. In C++, the type system is also a large part of the interface design.

For scientific code, types may encode distinctions such as:

```text
mesh
finite element
state vector
residual vector
control vector
quadrature rule
linear operator
solver configuration
```

A type distinction can prevent operations that would otherwise only fail numerically at run time. A strong type system does not prove the mathematics correct, but it can make many invalid program states harder to express.

### 1.2 `auto` asks the compiler to infer a type

C++ can often infer a type from an initializer:

```cpp
auto n = values.size();
auto layout = make_layout(dimensions);
```

`auto` does not make C++ dynamically typed. The compiler still determines one fixed static type for each variable.

It is most useful when the initializer already communicates the type clearly or when the exact type is long and mechanical. It can be harmful if the type itself carries information the reader needs. When reading code, mentally substitute “the compiler-deduced static type” rather than “some arbitrary type.”

## 2. Source files, headers, namespaces, and scope

C++ programs are usually divided across many files. The most common distinction is between:

- **headers**, often ending in `.h`, `.hpp`, or similar, which expose declarations and reusable definitions;
- **source files**, usually `.cc` or `.cpp`, which contain implementation code compiled into object files.

A header may contain an interface such as

```cpp
class LinearSolver
{
public:
  void solve(Vector &x, const Vector &b) const;
};
```

while a source file supplies the function body.

Templates are an important exception: their definitions usually need to be visible where the compiler instantiates them, so template-heavy scientific libraries often place substantial implementation code in headers.

The build consequences of headers and source files are developed in [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md). For now, the important reading rule is that an `#include` makes declarations from another header available in the current translation unit.

### 2.1 Namespaces organize names

A library can place its names inside a namespace:

```cpp
namespace solver
{
  class Control;
}
```

The fully qualified name is

```cpp
solver::Control
```

The `::` operator means “look for this name inside that scope.” Nested namespaces are common:

```cpp
nmopt::contract::ExecutableModel
```

A namespace is not an object and does not have a lifetime. It is a naming mechanism.

### 2.2 Scope controls where a name is visible

Braces usually introduce a scope:

```cpp
if (converged)
  {
    const double residual = compute_residual();
    // residual is visible here
  }
// residual is no longer visible here
```

Scope and lifetime often align for ordinary local objects: a local object is normally destroyed when execution leaves its scope. This fact is the basis of RAII later in the chapter.

## 3. Values, references, and pointers answer different questions

Suppose we have

```cpp
Vector x;
```

Then `x` is an object with its own lifetime. There are several ways to give another function access to it.

### 3.1 Passing by value creates a separate parameter object

```cpp
void normalize(Vector x);
```

The parameter `x` is a separate object initialized from the caller's argument. For a copyable lvalue, this ordinarily means a copy. Changes to the local `x` do not modify the caller's vector.

Passing by value is natural for small scalar types:

```cpp
void set_tolerance(double tolerance);
```

It is also useful when the function intentionally wants its own object, possibly by moving from the caller.

### 3.2 A reference is another name for an existing object

```cpp
void normalize(Vector &x);
```

The `&` here declares an lvalue reference. The parameter refers to the caller's object rather than creating a second vector. Changes through `x` change that original object.

A reference must be bound to an object when it is created and cannot later be reseated to refer to a different object.

For read-only access, use a reference to const:

```cpp
double norm(const Vector &x);
```

This avoids copying the vector while promising not to modify it through `x`.

A useful first approximation for scientific interfaces is:

```text
small value                         → pass by value
large read-only object              → const reference
object intentionally modified       → non-const reference
```

This is a reading heuristic, not a complete function-design rule.

### 3.3 A raw pointer can express nullable access or low-level relationships

```cpp
void inspect(const Vector *x);
```

A pointer stores an address. Unlike a reference, it can be null:

```cpp
if (x != nullptr)
  use(*x);
```

`*x` dereferences the pointer to access the pointed-to object. `x->size()` is shorthand for `(*x).size()` when accessing a member.

A raw pointer does **not**, by itself, tell the reader who owns the object or how long it will live. Modern C++ code therefore tends to use raw pointers mainly for non-owning or low-level relationships and smart pointers when ownership must be represented explicitly.

### 3.4 Reference versus pointer is not just syntax

These two interfaces communicate different expectations:

```cpp
void apply(const Operator &A, Vector &y);
void apply(const Operator *A, Vector &y);
```

The first says that an operator object is required. The second allows at least the *representation* of “no operator” through `nullptr`, whether or not that is semantically accepted by the function.

When reading scientific code, ask first whether the relationship is required and non-owning, optional and non-owning, or ownership-carrying. The syntax often reflects that distinction.

## 4. `const` is an interface promise about mutation

The most common use of `const` in scientific code is the read-only reference:

```cpp
void assemble(const Mesh &mesh);
```

Inside `assemble`, the program cannot modify `mesh` through this reference, except through deliberately `mutable` internals that are outside our scope here.

This matters for reasoning: a caller can pass a large object without copying it while knowing the called function cannot directly change its logical state through that parameter.

### 4.1 A const member function promises not to modify the logical object through `this`

Consider

```cpp
class VectorSpace
{
public:
  std::size_t dimension() const;
};
```

The trailing `const` belongs to the member function. It means that calling `dimension()` does not modify the object through the ordinary non-mutable data members.

Therefore it can be called on a const object:

```cpp
const VectorSpace space;
auto n = space.dimension();
```

This is different from a return type being const.

### 4.2 Const pointers have several forms

These declarations are worth recognizing even if you rarely write all of them:

```cpp
const Vector *p1;       // pointer to const Vector
Vector *const p2 = ...; // const pointer to mutable Vector
const Vector *const p3 = ...; // const pointer to const Vector
```

Read from the identifier outward:

- `p1` may point somewhere else later, but cannot modify the `Vector` through `p1`;
- `p2` always stores the same address, but may modify that `Vector`;
- `p3` can do neither.

In higher-level scientific interfaces, references and smart pointers usually make ownership/lifetime intentions clearer than const raw-pointer combinations.

## 5. Classes collect state and operations

A class defines a user-defined type. A small example is

```cpp
class LinearModel
{
public:
  LinearModel(double coefficient)
    : coefficient_(coefficient)
  {}

  double
  evaluate(double x) const
  {
    return coefficient_ * x;
  }

private:
  double coefficient_;
};
```

The object stores `coefficient_` and exposes an operation `evaluate`.

`public` members form the visible interface. `private` members are implementation details accessible to the class itself and its friends.

A `struct` is almost the same language mechanism as a `class`; the main default difference is access control. By convention, a simple aggregate-like data record is often a `struct`, while a type maintaining invariants behind an interface is often a `class`.

### 5.1 Constructors establish usable objects

The declaration

```cpp
LinearModel(double coefficient)
```

is a constructor. It creates a `LinearModel` object from a coefficient.

The syntax

```cpp
: coefficient_(coefficient)
```

is a **member-initializer list**. It initializes the data member before the constructor body executes.

This is not equivalent to first default-constructing every member and then assigning to it. Some members, such as references and const objects, must be initialized directly, and many resource-owning types are more naturally constructed this way.

### 5.2 Invariants explain why constructors matter

Suppose a matrix-free operator is only meaningful when its mesh and finite-element description are compatible. A constructor can require those dependencies and reject invalid inputs. After construction, member functions may then assume the invariant holds.

This is a common scientific-library pattern:

```text
constructor
    establish required relationships

public operations
    rely on those relationships
```

A well-designed object therefore communicates more than a bag of values; it can represent a valid state of a numerical concept.

## 6. Lifetime: every reference must refer to an object that still exists

Consider

```cpp
const Vector &
make_bad_vector()
{
  Vector temporary;
  return temporary;
}
```

This is wrong. `temporary` is destroyed when the function returns, so the returned reference dangles.

The central lifetime rule is simple to state:

> A pointer, reference, view, or callback capture must never be used after the object it refers to has been destroyed.

The difficulty is that large scientific programs create long chains of objects that refer to other objects. Understanding lifetime is therefore more important than memorizing syntax.

### 6.1 Automatic lifetime follows scope

A local object such as

```cpp
Vector x;
```

is automatically destroyed at the end of its scope. Its destructor runs even when the scope is exited by an exception.

This deterministic destruction is one of the defining features of C++ resource management.

### 6.2 Dynamic lifetime outlives the creating scope when ownership is retained

An object allocated dynamically can outlive the scope that initiated its creation. Modern C++ normally represents this with an owning object such as `std::unique_ptr` or `std::shared_ptr`, rather than with manual `new` and `delete` scattered through application code.

We will introduce these after the RAII principle that motivates them.

## 7. RAII ties resource lifetime to object lifetime

RAII stands for **Resource Acquisition Is Initialization**. The name is historical; the practical idea is:

> Put the responsibility for releasing a resource inside an object's destructor, then let ordinary object lifetime control cleanup.

A resource may be:

- dynamically allocated memory;
- a file handle;
- a mutex lock;
- a temporary operating-system resource;
- a library object requiring a paired cleanup operation.

For example, `std::vector<double>` owns dynamic memory internally. You normally write

```cpp
{
  std::vector<double> values(1000);
  // use values
}
```

and do not manually free the memory. When `values` is destroyed, its destructor releases the owned storage.

### 7.1 RAII is broader than smart pointers

A common misconception is that RAII means “use smart pointers.” Smart pointers are one RAII tool. `std::vector`, `std::string`, file streams, locks, and many library classes also use the same principle.

The desired ownership graph is usually built from ordinary values first and explicit smart pointers only where dynamically shared or transferred ownership is actually required.

### 7.2 Why RAII matters in numerical code

Scientific programs can fail midway through setup because a mesh file is invalid, a matrix factorization fails, or a run-time contract is violated. With RAII, objects already constructed on the stack are still destroyed as the exception unwinds through scopes.

This makes resource management local and compositional rather than requiring every exit path to remember every cleanup action.

## 8. Ownership and borrowing are different relationships

Suppose an object `Solver` needs access to an `Operator`.

One possibility is that `Solver` merely borrows an operator managed elsewhere:

```cpp
class Solver
{
public:
  Solver(const Operator &op)
    : op_(op)
  {}

private:
  const Operator &op_;
};
```

Then `Solver` does not own the operator. The surrounding program must ensure that `op_` outlives the `Solver` use that depends on it.

Another possibility is unique ownership:

```cpp
class Solver
{
private:
  std::unique_ptr<Operator> op_;
};
```

Now the `Solver` owns an operator dynamically. Destroying the `Solver` destroys its owned operator.

A third possibility is genuinely shared ownership:

```cpp
std::shared_ptr<const Operator> op_;
```

Several objects may now share responsibility for keeping the same immutable operator alive.

These three designs have different semantics. The choice should follow the ownership relationship rather than convenience.

## 9. `std::unique_ptr`: one owner, transferable ownership

A `std::unique_ptr<T>` represents exclusive ownership of a dynamically allocated `T`.

```cpp
auto matrix = std::make_unique<Matrix>();
```

Only one `unique_ptr` owns that object at a time. Copying the pointer is forbidden because that would create two exclusive owners. Ownership can instead be transferred:

```cpp
auto other = std::move(matrix);
```

After the transfer, `other` owns the object. The moved-from `matrix` no longer owns that resource.

### 9.1 Exclusive ownership is often the simplest ownership

If one subsystem creates a resource and no other subsystem needs to extend its lifetime independently, `unique_ptr` communicates this clearly.

A function can transfer ownership into another object:

```cpp
void set_preconditioner(std::unique_ptr<Preconditioner> p);
```

The caller must then move its pointer:

```cpp
solver.set_preconditioner(std::move(p));
```

The syntax makes the ownership transfer visible.

## 10. `std::shared_ptr`: shared responsibility for lifetime

A `std::shared_ptr<T>` keeps an object alive while at least one owning shared pointer still refers to it.

```cpp
auto layout = std::make_shared<BlockLayout>(...);
```

Copying the `shared_ptr` creates another owner:

```cpp
auto another_view = layout;
```

Both now share responsibility for the lifetime of the same `BlockLayout` object.

### 10.1 Shared ownership is not “a safer pointer” by default

It solves a specific lifetime problem: several components genuinely need to retain ownership independently.

Using `shared_ptr` everywhere can obscure who logically owns what and introduces reference-counting machinery. Prefer ordinary values or unique ownership when the ownership graph is simpler.

### 10.2 `shared_ptr<const T>` is a useful scientific-interface pattern

Suppose multiple vectors or operators must refer to the same structural layout, but none should modify it. Then

```cpp
std::shared_ptr<const BlockLayout>
```

expresses two facts:

- the layout has shared lifetime;
- users access it as immutable.

This is useful when many numerical values must remain consistent with one structural object.

### 10.3 Cycles motivate `std::weak_ptr`

If object A owns B through a `shared_ptr` and B owns A through another `shared_ptr`, their reference counts may never reach zero.

`std::weak_ptr<T>` observes an object managed by `shared_ptr` without contributing to ownership. It is useful for breaking such cycles or representing optional non-owning access into a shared-lifetime graph.

For the background needed here, it is enough to recognize this role; we do not need the full weak-pointer API.

## 11. Move semantics transfer resources without requiring a deep copy

A large matrix, vector, or callback object may own resources that are expensive to duplicate. C++ distinguishes copying from moving.

Suppose a constructor receives an object by value:

```cpp
class Model
{
public:
  Model(std::string label)
    : label_(std::move(label))
  {}

private:
  std::string label_;
};
```

The parameter `label` is a local object. Once its contents have been transferred into `label_`, the constructor no longer needs the old value. `std::move(label)` says that the member initialization may use move operations rather than requiring a copy.

### 11.1 `std::move` does not itself move resources

This distinction is important. `std::move(x)` is essentially a cast that marks `x` as eligible to bind to move-aware operations. The actual transfer is performed by the move constructor or move assignment operator selected by overload resolution.

So this:

```cpp
auto y = std::move(x);
```

means roughly:

```text
allow y's construction to treat x as a source whose resources may be transferred
```

not:

```text
perform a universal magic move operation on x
```

### 11.2 What can be done with a moved-from object?

For standard-library objects, a moved-from value is generally left in a valid but unspecified state unless a more specific guarantee is documented. It can be destroyed or assigned a new value, and operations whose preconditions are known to hold remain valid. You should not assume that it still contains its old value.

For a user-defined type, the type's documented move contract matters.

A useful reading habit is: after `std::move(x)`, look for whether `x` is used again. If it is, verify that the use is valid rather than assuming the original content survived.

## 12. Templates parameterize code by types

Scientific libraries often need the same algorithm to work with several vector, matrix, scalar, or dimension types. A class template expresses this at compile time:

```cpp
template <typename Vector>
class LinearOperator
{
public:
  void apply(const Vector &x, Vector &y) const;
};
```

`Vector` is a template parameter. `LinearOperator<MyVector>` and `LinearOperator<OtherVector>` are distinct C++ types generated from the same template definition.

### 12.1 Templates are compile-time abstraction

A template does not mean that a run-time variable stores a type. The compiler instantiates code for the concrete template arguments that appear in the program.

This is why scientific C++ can write generic code without necessarily paying for a run-time virtual dispatch on every operation.

A common deal.II-style example is a dimension template:

```cpp
template <int dim>
class Problem;
```

Then `Problem<2>` and `Problem<3>` are different types. The dimension can influence array sizes and select dimension-specific code at compile time.

### 12.2 Type aliases shorten long names

Inside a template one often finds

```cpp
using Vector = typename Backend::Vector;
```

This gives a shorter local name to a type associated with `Backend`.

The keyword `typename` tells the compiler that the dependent name `Backend::Vector` is a type. Because `Backend` is a template parameter, the compiler cannot know that merely from parsing the template definition.

You do not need to master all dependent-name rules to read `nmopt`. The useful interpretation is simply:

```text
Backend supplies an associated type called Vector
```

### 12.3 Template errors can be long because the compiler reports the instantiation chain

If a generic function is instantiated with a type that lacks an expected operation, the error may include several nested template frames. When reading such diagnostics, find:

1. the first project-owned line where the template was instantiated;
2. the concrete template arguments;
3. the innermost meaningful complaint, such as “no matching function” or “no member named ...”.

The build-workflow companion develops diagnostic reading more generally.

## 13. Virtual interfaces provide run-time polymorphism

Templates choose a type at compile time. Sometimes a program needs to treat several implementations through one run-time interface instead.

Consider

```cpp
class Model
{
public:
  virtual ~Model() = default;

  virtual double objective(const Vector &x) const = 0;
  virtual Vector residual(const Vector &x) const = 0;
};
```

The `= 0` makes these functions **pure virtual**. `Model` is therefore abstract: it defines an interface but cannot itself be instantiated as a complete implementation.

A derived class supplies the operations:

```cpp
class PoissonModel final : public Model
{
public:
  double objective(const Vector &x) const override;
  Vector residual(const Vector &x) const override;
};
```

### 13.1 `override` lets the compiler check the intended relationship

If the signature does not actually override a virtual base function, `override` produces a compile-time error. This catches mistakes such as changing a `const` qualifier or parameter type accidentally.

### 13.2 `final` closes further overriding or inheritance

On a class,

```cpp
class PoissonModel final : public Model
```

means that no class may derive from `PoissonModel`.

On a virtual function, `final` means that later derived classes may not override that function.

### 13.3 Why an interface base class normally has a virtual destructor

If a derived object is owned through a base pointer and deleted through that base interface, destruction must dispatch correctly to the derived destructor. A virtual destructor provides that behavior.

This is why interface classes commonly begin with

```cpp
virtual ~Model() = default;
```

### 13.4 Templates and virtual interfaces solve different problems

A rough comparison is:

```text
template parameter
    implementation selected at compile time
    concrete type appears in the static type

virtual interface
    implementation selected through a run-time object
    callers can use a common base interface
```

Large scientific libraries often combine both. For example, an abstract interface itself may be templated on the vector backend.

## 14. Lambdas package a small piece of behavior as an object

A lambda expression creates an unnamed callable object:

```cpp
auto square = [](double x) {
  return x * x;
};

const double y = square(3.0);
```

The empty `[]` is the **capture list**. It tells the compiler which variables from the surrounding scope are stored or referenced by the lambda object.

### 14.1 Capture by value stores a copy in the lambda object

```cpp
double alpha = 2.0;
auto scale = [alpha](double x) {
  return alpha * x;
};
```

The lambda contains its own captured value of `alpha`. Changing the outer variable later does not change that stored value.

### 14.2 Capture by reference borrows the surrounding object

```cpp
double alpha = 2.0;
auto scale = [&alpha](double x) {
  return alpha * x;
};
```

Now the lambda refers to the original `alpha`. If that lambda is stored and invoked after `alpha` has been destroyed, the capture dangles.

This is the same lifetime issue as any other reference.

### 14.3 `[=]` and `[&]` are broad capture defaults

```cpp
[=](double x) { ... } // implicitly capture needed surrounding variables by value
[&](double x) { ... } // implicitly capture needed surrounding variables by reference
```

They are concise but can hide what a callback retains. In long-lived scientific callbacks, explicit captures often make lifetime reasoning easier.

### 14.4 Capturing `this` refers to the surrounding object

Inside a member function, a lambda may access members of the current object. In C++17, capturing `this` gives the lambda access to that object rather than making an independent deep copy of all members.

Therefore a stored callback that refers to `this` must not outlive the object whose members it uses.

This matters when an application registers callbacks into a solver or compiler object that may retain them.

## 15. `std::function` erases the concrete callable type

Each lambda expression has its own compiler-generated type. A function pointer has another type. A function object class has another. Sometimes an interface wants to accept any of these as long as they share the same call signature.

C++17 provides `std::function` for this purpose:

```cpp
using ResidualAction =
  std::function<Vector(const Vector &)>;
```

A `ResidualAction` can hold many different callable targets:

```cpp
Vector residual_function(const Vector &x);

ResidualAction a = residual_function;
ResidualAction b = [](const Vector &x) {
  return compute_residual(x);
};
```

Both can later be invoked uniformly:

```cpp
Vector r = a(x);
```

### 15.1 This is a form of type erasure

The caller knows the required interface

```text
Vector(const Vector &)
```

without knowing the exact concrete type of the stored callable.

This is useful at application boundaries. A generic solver does not need to know whether a residual comes from:

- a free function;
- a lambda capturing an application object;
- a small function-object class;
- a wrapper around a library call.

It only needs the callable contract.

### 15.2 Type erasure has a cost, but do not guess the implementation cost

Compared with an inlined concrete lambda type, `std::function` adds an indirection layer and stores a polymorphic callable. Implementations may use small-object optimizations, and some targets may require dynamic allocation. The exact cost is implementation- and target-dependent.

For interfaces where one callback invocation triggers a PDE residual assembly or solve, this wrapper overhead is often negligible relative to the numerical work. In a tiny inner loop called billions of times, the tradeoff would deserve separate analysis.

### 15.3 An empty `std::function` is a distinct state

A default-constructed `std::function` contains no callable target. It can be checked in a Boolean context:

```cpp
if (residual)
  residual(x);
```

Invoking an empty `std::function` throws `std::bad_function_call`. Scientific APIs often validate required callbacks during construction so that missing operations fail at the boundary rather than much later in an algorithm.

## 16. A callback-based scientific interface combines the previous ideas

Consider a minimal generic model:

```cpp
#include <functional>
#include <utility>


template <typename Vector>
class CallbackModel
{
public:
  using Residual = std::function<Vector(const Vector &)>;
  using LinearAction =
    std::function<Vector(const Vector &, const Vector &)>;

  CallbackModel(Residual residual,
                LinearAction jvp,
                LinearAction vjp)
    : residual_(std::move(residual))
    , jvp_(std::move(jvp))
    , vjp_(std::move(vjp))
  {}

  Vector
  residual(const Vector &x) const
  {
    return residual_(x);
  }

  Vector
  jvp(const Vector &x, const Vector &direction) const
  {
    return jvp_(x, direction);
  }

  Vector
  vjp(const Vector &x, const Vector &seed) const
  {
    return vjp_(x, seed);
  }

private:
  Residual residual_;
  LinearAction jvp_;
  LinearAction vjp_;
};
```

There are several layers of meaning here.

`Vector` is a template parameter, so the model is not tied to one concrete vector backend.

The aliases `Residual` and `LinearAction` assign readable names to callable signatures.

The constructor receives callback objects by value. It then moves them into the stored members. This is a common “take ownership of the supplied value” idiom.

The member functions accept their inputs as const references, so large vector objects need not be copied just to evaluate an operation.

The callbacks are stored behind `std::function`, so the caller can supply different concrete implementations without deriving a new class from `CallbackModel`.

### 16.1 Supplying callbacks with lambdas

Suppose a matrix-like object `A` provides `vmult` and `Tvmult` operations. An application could write schematically:

```cpp
CallbackModel<Vector> model(
  [&A](const Vector &x) {
    Vector r(x.size());
    A.vmult(r, x);
    return r;
  },
  [&A](const Vector &, const Vector &d) {
    Vector result(d.size());
    A.vmult(result, d);
    return result;
  },
  [&A](const Vector &, const Vector &q) {
    Vector result(q.size());
    A.Tvmult(result, q);
    return result;
  });
```

The lambda callbacks borrow `A` by reference. Therefore `A` must remain alive for as long as `model` might invoke those callbacks.

Changing the capture to a shared owning pointer could change the lifetime relationship:

```cpp
auto A = std::make_shared<Operator>(...);

auto apply = [A](const Vector &x) {
  // the lambda now retains shared ownership of A
};
```

This is a design decision, not merely syntactic variation.

## 17. `enum class` gives named alternatives without an open-ended integer convention

Suppose a solver can choose one of several strategies:

```cpp
enum class SolverKind
{
  cg,
  gmres,
  minres
};
```

Then an interface can accept

```cpp
SolverKind method
```

rather than a bare integer such as `0`, `1`, or `2`.

Scoped enumerations avoid implicitly polluting the surrounding namespace with their enumerator names and do not implicitly convert to integers in the same loose way as old unscoped enums.

They are especially useful for small, closed sets of semantic choices.

## 18. `std::optional` represents “a value may be absent” explicitly

An optional configuration can be modeled as

```cpp
std::optional<double> absolute_tolerance;
```

The object either contains a `double` or contains no value.

Typical operations include

```cpp
if (absolute_tolerance)
  use(*absolute_tolerance);
```

or

```cpp
const double tol = absolute_tolerance.value_or(default_tolerance);
```

This is often clearer than inventing a sentinel such as `-1.0` to mean “not supplied,” especially if negative numbers could otherwise be meaningful or invalid for an unrelated reason.

For our purposes, recognize `optional<T>` as **zero or one `T`**, not as a pointer ownership mechanism.

## 19. `std::variant` represents one value chosen from a closed set of types

A variant can store exactly one of several listed types:

```cpp
using BoundaryData =
  std::variant<ConstantData, FunctionData, TableData>;
```

Unlike a base-class pointer, this set of alternatives is closed in the variant's type itself. The program can inspect or visit the active alternative.

Scientific configuration and semantic models sometimes use variants when a field can have one of several structurally different representations.

A full treatment of `std::visit`, exhaustive visitors, and algebraic-data-type style programming is beyond our needs. The main reading rule is:

```text
variant<A, B, C>
    exactly one currently active value of type A, B, or C
```

## 20. Exceptions and contracts communicate failure differently from return values

A scientific interface often has preconditions. For example, two block vectors may need compatible layouts, or a requested solver may require an SPD operator.

One possible boundary check is

```cpp
if (!compatible(left, right))
  throw std::invalid_argument("incompatible layouts");
```

An exception transfers control to an enclosing handler rather than returning a normal value.

### 20.1 Exceptions do not replace numerical convergence reports

There is an important distinction between:

```text
programming / contract failure
    invalid argument
    impossible configuration
    missing required callback

algorithmic outcome
    nonlinear solver did not converge
    maximum iterations reached
    tolerance not met
```

A numerical library may choose to represent the second category with a structured report rather than with an exception, because nonconvergence can be an expected outcome that the caller needs to inspect.

The exact policy is library-specific. When reading code, determine whether an exception indicates a broken interface contract, an exceptional run-time condition, or an algorithmic failure translated into exception form.

### 20.2 RAII makes exceptions manageable

Because local resource-owning objects are destroyed during stack unwinding, throwing an exception does not require handwritten cleanup for every ordinary RAII-managed resource.

This connection between error handling and lifetime is one reason RAII is central to modern C++ rather than an isolated smart-pointer technique.

## 21. Reading scientific C++: separate language mechanics from numerical meaning

Consider

```cpp
template <typename Backend>
CovectorBlockT<Backend>
residual_vjp(const PrimalBlockT<Backend> &x,
             const PrimalBlockT<Backend> &seed) const;
```

The language mechanics say:

- this is part of a template depending on `Backend`;
- it returns a covector-block object specialized for that backend;
- both arguments are borrowed read-only;
- the member function is const.

The numerical meaning may be:

- `x` is the point where the residual is linearized;
- `seed` is a primal representation in the residual test space;
- the return value represents the transpose/adjoint action in the variable dual space.

Neither level replaces the other. A reader needs enough C++ to get past the mechanics, then should return attention to the numerical contract.

This separation is especially useful when code contains long type names. First decode what the type system promises; then ask what mathematical object each role represents.

## 22. Where these patterns appear in `nmopt`

The current `nmopt` public contract is a useful example because it combines most of the mechanisms in this chapter without requiring advanced metaprogramming.

The project declares its public interface as C++17 in the root CMake configuration. Its callback executable model is a class template parameterized by a backend. It defines aliases such as conceptually

```cpp
using ResidualAction =
  std::function<Covector(const Primal &)>;
```

and stores residual, JVP, VJP, objective, and objective-derivative callbacks. The constructor receives these values and moves them into members.

The same type derives from an abstract executable-model interface whose operations are pure virtual. Implementations mark their methods with `override`, and callback models are `final` because they are intended as complete concrete adapters rather than intermediate base classes.

The block-layout layer uses

```cpp
std::shared_ptr<const BlockLayout>
```

to share immutable structural information among block values. Its generic block types use a backend-associated vector type through declarations of the form

```cpp
typename Backend::Vector
```

These are direct instances of the lifetime, ownership, templates, aliases, and run-time interface ideas developed above.

The project-specific semantics of those classes belong to the `nmopt` manual. The point of this chapter is that their C++ form should no longer obscure the numerical roles.

## 23. Scope frontier

This chapter deliberately omits large parts of C++. They are not needed to follow the present `nmopt`/deal.II interface path.

We stop before:

- template metaprogramming and SFINAE in detail;
- C++20 concepts and ranges;
- custom allocators and memory-resource APIs;
- lock-free programming, threads, and concurrency primitives;
- low-level placement construction and manual lifetime management;
- complex overload-resolution rules;
- perfect forwarding and forwarding-reference design;
- advanced iterator and ranges programming;
- custom exception hierarchies and formal exception-safety guarantees;
- ABI details and dynamic libraries;
- compiler optimization internals;
- C++ module systems.

You may encounter some of these in larger C++ ecosystems. When that happens, learn the concept because a concrete dependency demands it rather than trying to master the language in advance.

## Where to go next

- **If build tools are unfamiliar:** [08a · Scientific software workflow: shell, compilation, CMake, and tests](08a-scientific-software-workflow.md) explains how source files become executables through compilation, linking, CMake, generated build systems, and tests.
- **If you already build C++ projects comfortably:** [09 · The deal.II finite-element workflow](09-dealii-finite-element-workflow.md) maps [03 · Finite elements](03-finite-elements.md) onto deal.II's object model using the C++ vocabulary developed here.
- **Into the `nmopt` manual:** the most relevant pages are:
  - [manual 11 · Compilation and lowering](../concepts/11-compilation-and-lowering.md), for callbacks, type-erased services, ownership, and compiler-produced numerical objects;
  - [manual 12 · Authoring and using compiled problems](../concepts/12-authoring-and-using-compiled-problems.md), for the application-facing compiled workflow;
  - [manual 13 · Integrating an existing PDE application](../concepts/13-integrating-an-existing-pde-application.md), for adapting application-owned objects and lifetimes to `nmopt` contracts.

## References and further reading

Publication details, stable links, and access notes are collected in the [Background reference catalogue](references.md).

The C++ suggestions below are optional and **have not been vetted as part of this background sequence**. They are possible continuations, not prerequisites.

- **Bjarne Stroustrup, *A Tour of C++*, third edition.**
  - Compact overview for readers who already have programming experience; substantially broader than the C++ subset needed here.
- **Bjarne Stroustrup and Herb Sutter, *C++ Core Guidelines*.**
  - Living design-guideline reference rather than a sequential tutorial.
  - The material on interfaces, parameter passing, classes, resource management, constants, templates, source files, and the standard library is especially relevant here.
- **cppreference.**
  - Precise lookup reference for language rules and facilities such as smart pointers, `std::move`, lambdas, `std::function`, `std::optional`, and `std::variant`.
  - Best used after a concept has been introduced rather than as a first textbook.
- **Wolfgang Bangerth, MATH 676 video lectures.**
  - Not a general C++ course, but the lectures on templates, scientific-program structure, debug versus optimized builds, and large scientific-software development provide useful context.
  - Older recordings can contain outdated deal.II APIs; use current deal.II documentation as the authority when interfaces differ.
