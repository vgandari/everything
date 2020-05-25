---
title: CS Terms
---

glob pattern
:    ??

helper class
:    used to assist in providing some functionality, which isn’t the
     main goal of the application or class in which it is used; utility
     class is a special case of a helper class in which the methods are
     all static https://en.wikipedia.org/wiki/Helper_class

hash table
:    a data structure that implements an associative array abstract data
     type, a structure that can map keys to values
     https://en.wikipedia.org/wiki/Hash_table

open addressing
:    https://en.wikipedia.org/wiki/Open_addressing 

linear probing
:    https://en.wikipedia.org/wiki/Linear_probing 

referent of the pointer 
:    ??
 
RAII
:    ??

smart pointer
:    abstract data type that simulates a pointer while providing added
     features, such as automatic memory management or bounds checking.
     [https://en.wikipedia.org/wiki/Smart_pointer](https://en.wikipedia.org/wiki/Smart_pointer)

reference counting
:    ??

circular references
:    ??

shared pointer
:    ??

memory pool
:    ??

type safety
:    ?? 

thread model
:    ?? e.g. POSIX

thread safety
:    Thread-safe code only manipulates shared data structures in a
     manner that ensures that all threads behave properly and fulfill
     their design specifications without unintended interaction.

race condition
:    ?? 

concurrency control
:    ?? 

critical section
:    section of a program where accessing a shared resource can lead to unexpected or erroneous behavior; cannot be executed by more than one process at a time. 

mutex (mutual exclusion)
:    requirement that one thread of execution never enters its critical
     section at the same time that another concurrent thread of
     execution enters its own critical section;
     https://en.wikipedia.org/wiki/Mutual_exclusion
     
closure
:    [https://en.wikipedia.org/wiki/Closure_(computer_programming)](https://en.wikipedia.org/wiki/Closure_(computer_programming))

## Python

[Common Object Structures](https://docs.python.org/3/c-api/structures.html)

dictionary
:    The goal of a hash function is to distribute the keys evenly in the
     array. A good hash function minimizes the number of collisions e.g.
     different keys having the same hash. Python does not have this kind
     of hash function. Python dictionaries use open addressing with
     linear probing to prevent collisions
     https://www.laurentluce.com/posts/python-dictionary-implementation/

map
:    ?? 

:::{.alert .alert-info}
**NOTE:**

It looks like Python dictionaries and C++ `std::unordered_map` are
implemented as hash tables
(https://stackoverflow.com/questions/8265608/why-is-pythons-dict-implemented-as-hash-table-whereas-stdmap-is-tree-based).
It also looks like `std::map` is slower.

:::

:::{.alert .alert-info}
**NOTE:**

A C++ `std::map` and a Python `map` are _not_ the same thing.

:::

global interpretter lock
:    prevents parallel execution https://wiki.python.org/moin/GlobalInterpreterLock, Understanding the Python GIL 

## Incomplete

associative array
hook
interrupt vectors
page fault
register
stack vs heap
dynamic, static, and automatic memory allocation
mask

lvalue
rvalue
object-oriented 
function pointers, also called ‘closures’ or ‘delegates’

language lawyer
checked iterator
in-cache/out-of-cache
register
[pipeline](https://en.wikipedia.org/wiki/Pipeline_(computing))
instruction pointer (IP)
program counter (PC)
RISC architecture
CISC architecture
branch instruction
out-of-order CPU
[branch prediction](https://en.wikipedia.org/wiki/Branch_predictor)
Vtable Run-Time Dispatch
## Notes

every list slicing operation involves making a copy of the relevant object references (but not the objects themselves).

