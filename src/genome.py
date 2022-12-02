"""A circular genome for simulating transposable elements."""

from typing import (
    Generic, TypeVar, Iterable,
    Callable, Protocol
)

T = TypeVar('T')

from abc import (
    # A tag that says that we can't use this class except by specialising it
    ABC,
    # A tag that says that this method must be implemented by a child class
    abstractmethod
)

from collections import namedtuple


class Genome(ABC):
    """Representation of a circular enome."""

    def __init__(self, n: int):
        """Create a genome of size n."""
        ...  # not implemented yet

    @abstractmethod
    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.

        Insert a new transposable element at position pos and len
        nucleotide forward.

        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.

        Returns a new ID for the transposable element.
        """
        ...  # not implemented yet

    @abstractmethod
    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.

        Copy the transposable element te to an offset from its current_node
        location.

        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.

        If te is not active, return None (and do not copy it).
        """
        ...  # not implemented yet

    @abstractmethod
    def disable_te(self, te: int) -> None:
        """
        Disable a TE.

        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ...  # not implemented yet

    @abstractmethod
    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        ...  # not implemented yet

    @abstractmethod
    def __len__(self) -> int:
        """Get the current_node length of the genome."""
        ...  # not implemented yet

    @abstractmethod
    def __str__(self) -> str:
        """
        Return a string representation of the genome.

        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.

        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        ...  # not implemented yet


class ListGenome(Genome):
    """
    Representation of a genome.
    Implements the Genome interface using Python's built-in lists
    """

    def __init__(self, n: int):
        """Create a new genome with length n."""
        #initialize the genome with no TE's yet.
        self.nucleotide= ['-'] * n
        self.TE = {}
        self.active_TE = []

    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.
        Insert a new transposable element at position pos and len
        nucleotide forward.
        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.
        Returns a new ID for the transposable element.
        """

        if len(self.TE) != 0:
            TE_ID = max(self.TE) + 1
        else:
            TE_ID = 1
        self.TE[TE_ID] = length
        self.active_TE.append(TE_ID)
        if self.nucleotide[pos] in self.active_TE:
            self.disable_te(self.nucleotide[pos])
        self.nucleotide[pos:pos] = [TE_ID]*length
        return TE_ID


    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.
        Copy the transposable element te to an offset from its current_node
        location.
        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.
        If te is not active, return None (and do not copy it).
        """
        if te not in self.active_TE:
            return None
        pos = self.nucleotide.index(te)
        res = (pos + offset) % len(self.nucleotide)
        return self.insert_te(res, self.TE[te])
        
            

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.
        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        self.active_TE.remove(te)
        for i in range(len(self)):
            if self.nucleotide[i] == te:
                for j in range(self.TE[te]):
                    self.nucleotide[i+j] = 'x'
                break

    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return self.active_TE

    def __len__(self) -> int:
        """current_node length of the genome."""
        return len(self.nucleotide)

    def __str__(self) -> str:
        """
        Return a string representation of the genome.
        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.
        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        nucleotide = ""
        for x in self.nucleotide:
            if isinstance(x, int):
                nucleotide += "A"
            else: 
                nucleotide += str(x)
        return nucleotide

class Link(Generic[T]):
    """Doubly linked link.""" 

    def __init__(self, val: T, p, n):
        """Create a new link and link up prev and next."""
        self.val = val
        self.prev = p
        self.next = n


def insert_next(link: Link[T], val: T) -> None:
    """Add a new link containing avl after link."""
    new_link = Link(val, link, link.next)
    new_link.prev.next = new_link
    new_link.next.prev = new_link

Ft = namedtuple("Ft", ["ft", "length"])


class LinkedListGenome(Genome):
    """
    Representation of a genome.
    Implements the Genome interface using linked lists.
    """
    head: Link[Ft]  
    active = {}
    counter = 0

    def __init__(self, n: int):
        """Create a new genome with length n."""
        self.head = Link(None, None, None)
        self.head.next = self.head
        self.head.prev = self.head
        insert_next(self.head.prev, Ft(0, n))

    def insert_te(self, pos: int, length: int) -> int:
        """
        Insert a new transposable element.
        Insert a new transposable element at position pos and len
        nucleotide forward.
        If the TE collides with an existing TE, i.e. genome[pos]
        already contains TEs, then that TE should be disabled and
        removed from the set of active TEs.
        Returns a new ID for the transposable element.
        """
        for ft, end in self.iter_pos():
            if end <= pos:
                continue
            else:
                half_1 = ft.val.length -end + pos
                half_2 = end - pos
                self.insert_te_at_pos(length, ft, half_1, half_2)
                return self.counter
        return None

    def insert_te_at_pos(self, length, ft, half_1, half_2):
        if  ft.val.ft == 1:
            self.disable_feature(ft)
        insert_next(ft, Ft(1, length))
        ft.val = Ft(ft.val.ft, half_1)
        insert_next(ft.next, Ft(ft.val.ft, half_2))
        self.counter += 1
        self.active[self.counter] = ft.next

    def disable_feature(self, ft):
        ft.val = Ft(2, ft.val.length)
        for identifier, disable in self.active.items():
            if disable is ft:
                self.active.pop(identifier)
                break


    def copy_te(self, te: int, offset: int) -> int | None:
        """
        Copy a transposable element.
        Copy the transposable element te to an offset from its current
        location.
        The offset can be positive or negative; if positive the te is copied
        upwards and if negative it is copied downwards. If the offset moves
        the copy left of index 0 or right of the largest index, it should
        wrap around, since the genome is circular.
        If te is not active, return None (and do not copy it).
        """

        ft = self.active[te]
        temp = self.active[te]
        
        if temp.val.ft != 1:
            return None

        if offset > 0:
            offset = offset - temp.val.length
            while offset > 0:
                ft = getattr(ft, "next")
                if ft is self.head:
                    ft = getattr(ft, "next")
                offset -= ft.val.length
            half_1 = offset + ft.val.length
            half_2 = abs(offset)
            self.insert_te_at_pos(temp.val.length, ft, half_1, half_2)
        else:
            offset *= -1
            while offset > 0:
                ft = getattr(ft, "prev")
                if ft is self.head:
                    ft = getattr(ft, "prev")
                offset -= ft.val.length
            self.insert_te_at_pos(temp.val.length, ft, abs(offset), offset + ft.val.length)
        return self.counter

    def disable_te(self, te: int) -> None:
        """
        Disable a TE.
        If te is an active TE, then make it inactive. Inactive
        TEs are already inactive, so there is no need to do anything
        for those.
        """
        ft = self.active.pop(te)
        if ft is not None:
            ft.val = Ft(2, ft.val.length)


    def active_tes(self) -> list[int]:
        """Get the active TE IDs."""
        return list(self.active.keys())

    def __len__(self) -> int:
        """Current length of the genome."""
        return len(self.iter_pos)

    def __str__(self) -> str:
        """
        Return a string representation of the genome.
        Create a string that represents the genome. By nature, it will be
        linear, but imagine that the last character is immidiatetly followed
        by the first.
        The genome should start at position 0. Locations with no TE should be
        represented with the character '-', active TEs with 'A', and disabled
        TEs with 'x'.
        """
        nucleotide = ""
        for element in self:
            if  element.ft == 0:
                nucleotide += element.length*"-"
            elif element.ft == 1:
                nucleotide += element.length*"A"
            else:
                nucleotide += element.length*"x"
        return nucleotide

    def iter_pos(self):
        count = 0
        ft = self.head.next
        while ft is not self.head:
            count += ft.val.length 
            yield (ft, count)
            ft = ft.next

    def __iter__(self):
        ft = self.head.next
        while ft is not self.head:
            yield ft.val
            ft = ft.next


