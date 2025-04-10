use core::fmt::Debug;
use core::marker::PhantomData;
use core::mem;
use core::ops::{Index, IndexMut};
use core::ptr::NonNull;
use std::alloc::{self, Layout, alloc, dealloc, realloc};
use std::iter::FusedIterator;

// todo
// - store seglen instead of log_seglen?
// - extended testing
// - add shrinking the array to size
// - efficient item drop
// - c interop
// - size checks/debug assertions
// - clone trait
// - from slice/box trait
// - intoiter
// - insert
// - remove
// - zero sized types
// - benchmark random operations against vec
// - usize as index?
// - store last element ptr and segment details for faster push/pop?
// - ExactSizeIterator, size_hint(), DoubleEndedIterator, IntoIterator, FromIterator

/// A non-contiguos extensible array that achieves asymptotically optimal
/// space overhead of O(√n) while maintaining comparable time efficiency*.
/// Based on the data structure in this paper: <https://jmmackenzie.io/pdf/mm22-adcs.pdf>
pub struct SqVec<T> {
    dope: NonNull<NonNull<T>>,
    dope_cap: u32,
    len: u32,
    alloc_seg_count: u32,
    marker: PhantomData<T>,
}

#[macro_export]
macro_rules! sqvec {
    ( $( $x:expr ),* ) => {
        {
            let mut temp_vec = SqVec::new();
            $(
                temp_vec.push($x);
            )*
            temp_vec
        }
    };
}

impl<T> SqVec<T> {
    /// Constructs a new, empty `SqVec<T>`.
    ///
    /// The sqvec will not allocate until elements are pushed onto it.
    pub fn new() -> Self {
        assert!(
            core::mem::size_of::<T>() != 0,
            "This data structure does not support ZST's"
        );
        Self {
            dope: NonNull::dangling(),
            dope_cap: 0,
            len: 0,
            alloc_seg_count: 0,
            marker: PhantomData,
        }
    }

    /// calculates the dope_cap from `self.alloc_seg_count`
    fn dope_cap(&self) -> u32 {
        todo!()
    }

    /// grows the dope vector to the correct size to fit the next set of segments
    fn grow_dope(&mut self) {
        let (new_cap, new_layout) = if self.dope_cap == 0 {
            (2, unsafe {
                Layout::array::<NonNull<T>>(2).unwrap_unchecked()
            })
        } else {
            let new_cap = self.dope_cap + Self::segs_of_len(Self::lsln(self.len.saturating_sub(1)));
            let new_layout =
                Layout::array::<NonNull<T>>(new_cap as usize).expect("Allocation too large");
            (new_cap, new_layout)
        };

        let handle = if self.dope_cap == 0 {
            unsafe { alloc(new_layout) }
        } else {
            // this layout cannot error because its size is based on an allocation that already exists
            let old_layout =
                unsafe { Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap_unchecked() };
            let old_handle = self.dope.as_ptr() as *mut u8;
            unsafe { realloc(old_handle, old_layout, new_layout.size()) }
        };

        self.dope = NonNull::new(handle as *mut NonNull<T>)
            .unwrap_or_else(|| alloc::handle_alloc_error(new_layout));
        self.dope_cap = new_cap;
    }

    /// Allocates a new segment of length `1 << self.log_seglen` at index `self.alloc_seg_count` in dope.
    fn alloc_seg(&mut self) {
        // debug_assert_eq!(1 << Self::lsln(self.len), Self::lslnv(self.len));
        let layout =
            Layout::array::<T>(Self::lsln(self.len) as usize).expect("Allocation too large");
        let handle = unsafe { alloc(layout) };
        let ptr =
            NonNull::new(handle as *mut T).unwrap_or_else(|| alloc::handle_alloc_error(layout));
        unsafe { self.dope.add(self.alloc_seg_count as usize).write(ptr) };
        self.alloc_seg_count += 1;
    }

    // fn lslnv(v: u32) -> u32 {
    //     if v == 0 {
    //         2
    //     } else {
    //         v.next_power_of_two()
    //     }
    // }

    /// Adds `val` to the end of the SqVec, increasing its length by 1.
    pub fn push(&mut self, val: T) {
        let (segnum, offset) = Self::mapping(self.len);
        if offset == 0 && segnum == self.alloc_seg_count {
            if segnum == self.dope_cap {
                self.grow_dope();
            }
            self.alloc_seg();
        }
        unsafe { self.item_ptr_inner(segnum, offset).write(val) };
        self.len += 1;
    }

    /// Removes the last value on the SqVec and returns it in an Option.
    /// If the SqVec is empty, Returns None
    pub fn pop(&mut self) -> Option<T> {
        if self.len == 0 {
            None
        } else {
            self.len -= 1;
            Some(unsafe { self.item_ptr(self.len).read() })
        }
    }

    pub fn first(&self) -> Option<&T> {
        if self.len == 0 {
            None
        } else {
            Some(unsafe { self.dope.read().as_ref() })
        }
    }

    pub fn first_mut(&mut self) -> Option<&mut T> {
        if self.len == 0 {
            None
        } else {
            Some(unsafe { self.dope.read().as_mut() })
        }
    }

    pub fn last(&self) -> Option<&T> {
        if self.len == 0 {
            None
        } else {
            Some(unsafe { self.item_ptr(self.len - 1).as_ref() })
        }
    }

    pub fn last_mut(&mut self) -> Option<&mut T> {
        if self.len == 0 {
            None
        } else {
            Some(unsafe { self.item_ptr(self.len - 1).as_mut() })
        }
    }

    /// swaps the elements at index's ix1 and ix2
    pub fn swap(&mut self, ix1: u32, ix2: u32) {
        assert!(ix1 < self.len && ix2 < self.len, "Index out of range");
        unsafe { self.swap_unchecked(ix1, ix2) };
    }

    /// swaps the elements at index's ix1 and ix2 without checking ix1 or ix2 in range
    pub unsafe fn swap_unchecked(&mut self, ix1: u32, ix2: u32) {
        unsafe { mem::swap(self.item_ptr(ix1).as_mut(), self.item_ptr(ix2).as_mut()) };
    }

    /// Returns the length of the SqVec (number of elements)
    #[inline]
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Returns true if the sqvec is empty and false otherwise
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// clears the sqvec dropping all present elements but
    /// maintaining the currently allocated space
    pub fn clear(&mut self) {
        if mem::needs_drop::<T>() {
            for item in self.iter_mut() {
                unsafe { core::ptr::drop_in_place(item) };
            }
        }
        self.len = 0;
    }

    /// Maps the given index to the segment number and index within the
    /// segment (offset) where the element at index `v` can be found.
    #[inline]
    fn mapping(v: u32) -> (u32, u32) {
        let b = if v == 0 {
            1
        } else {
            (33 - v.leading_zeros()) >> 1
        };
        let segnum = (v >> b) + (1 << (b - 1)) - 1;
        let offset = v & ((1 << b) - 1);
        (segnum, offset)
    }

    /// Gives the dope index of the first segment of length seg_len
    #[inline]
    fn first_seg(seg_len: u32) -> u32 {
        (3 * seg_len >> 2) - 1
    }

    /// Returns the dope index of the last segment of length `seg_len`. `seg_len` must be a power of 2
    #[inline]
    fn last_seg(seg_len: u32) -> u32 {
        (3 * seg_len >> 1) - 2
    }

    /// Returns the max number of segments of length `seg_len`. `seg_len` must be a power of 2
    #[inline]
    fn segs_of_len(seg_len: u32) -> u32 {
        if seg_len == 2 { 2 } else { 3 * seg_len >> 2 }
    }

    /// Returns the pointer to the item in the SqVec at index `ix` (Does not check if ix is < self.len)
    #[inline]
    unsafe fn item_ptr(&self, ix: u32) -> NonNull<T> {
        let (segnum, offset) = Self::mapping(ix);
        unsafe { self.item_ptr_inner(segnum, offset) }
    }

    #[inline]
    unsafe fn item_ptr_inner(&self, segnum: u32, offset: u32) -> NonNull<T> {
        unsafe { self.dope.add(segnum as usize).read().add(offset as usize) }
    }

    /// Returns the segment length of the segment that would contain item of index `v`
    #[inline]
    fn lsln(v: u32) -> u32 {
        if v == 0 {
            2
        } else {
            1 << ((32 - v.leading_zeros()).div_ceil(2))
        }
    }

    fn raw_dope_iter(&self) -> RawDopeIter<T> {
        let (lseg, loff) = Self::mapping(self.len.saturating_sub(1));
        RawDopeIter {
            seg_ptr: if self.len == 0 {
                unsafe {
                    NonNull::new_unchecked(core::ptr::without_provenance_mut(
                        self.dope.addr().get() + 1,
                    ))
                }
            } else {
                self.dope
            },
            last_seg_ptr: unsafe { self.dope.add(lseg as usize) },
            last_offset: loff,
            seg_len: 2,
            segs_left: 2,
        }
    }

    fn raw_iter(&self) -> RawIter<T> {
        RawIter {
            rdi: self.raw_dope_iter(),
            rsi: {
                let ptr = NonNull::dangling();
                RawSegIter {
                    ptr: unsafe {
                        NonNull::new_unchecked(core::ptr::without_provenance_mut(
                            ptr.addr().get() + 1,
                        ))
                    },
                    end: ptr,
                }
            },
        }
    }

    /// Returns an iterator that returns immutable (shared) references to the
    /// SqVec's items in ascending order from 0..self.len
    pub fn iter(&self) -> Iter<T> {
        Iter {
            raw_iter: self.raw_iter(),
            marker: PhantomData,
        }
    }

    /// Returns an iterator that returns mutable references to the
    /// SqVec's items in ascending order from 0..self.len
    pub fn iter_mut(&mut self) -> IterMut<T> {
        IterMut {
            raw_iter: self.raw_iter(),
            marker: PhantomData,
        }
    }
}

struct RawDopeIter<T> {
    seg_ptr: NonNull<NonNull<T>>,
    last_seg_ptr: NonNull<NonNull<T>>,
    last_offset: u32,
    seg_len: u32,
    segs_left: u32,
}

struct RawSegIter<T> {
    ptr: NonNull<T>,
    end: NonNull<T>,
}

struct RawIter<T> {
    rdi: RawDopeIter<T>,
    rsi: RawSegIter<T>,
}

// is T: 'a necessary?
pub struct Iter<'a, T: 'a> {
    raw_iter: RawIter<T>,
    marker: PhantomData<&'a T>,
}

pub struct IterMut<'a, T: 'a> {
    raw_iter: RawIter<T>,
    marker: PhantomData<&'a mut T>,
}

impl<T> Iterator for RawDopeIter<T> {
    type Item = RawSegIter<T>;

    fn next(&mut self) -> Option<Self::Item> {
        let out = if self.seg_ptr < self.last_seg_ptr {
            let first_item = unsafe { self.seg_ptr.read() };
            Some(RawSegIter {
                ptr: first_item,
                end: unsafe { first_item.add(self.seg_len as usize - 1) },
            })
        } else if self.seg_ptr == self.last_seg_ptr {
            let first_item = unsafe { self.seg_ptr.read() };
            Some(RawSegIter {
                ptr: first_item,
                end: unsafe { first_item.add(self.last_offset as usize) },
            })
        } else {
            return None;
        };

        self.segs_left -= 1;
        if self.segs_left == 0 {
            self.seg_len <<= 1;
            self.segs_left = 3 * self.seg_len >> 2;
        }

        self.seg_ptr = unsafe { self.seg_ptr.add(1) };
        out
    }
}

impl<T> Iterator for RawSegIter<T> {
    type Item = NonNull<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr > self.end {
            return None;
        }

        let out = self.ptr;
        self.ptr = unsafe { self.ptr.add(1) };
        Some(out)
    }
}

impl<T> Iterator for RawIter<T> {
    type Item = NonNull<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.rsi.next() {
            Some(x)
        } else {
            if let Some(y) = self.rdi.next() {
                self.rsi = y;
                self.rsi.next()
            } else {
                None
            }
        }
    }
}

impl<'a, T: 'a> Iterator for Iter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_iter.next().map(|x| unsafe { x.as_ref() })
    }
}

impl<'a, T: 'a> Iterator for IterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        self.raw_iter.next().map(|mut x| unsafe { x.as_mut() })
    }
}

impl<T> FusedIterator for RawDopeIter<T> {}
impl<T> FusedIterator for RawSegIter<T> {}
impl<T> FusedIterator for RawIter<T> {}
impl<T> FusedIterator for Iter<'_, T> {}
impl<T> FusedIterator for IterMut<'_, T> {}

// Control flow could use some work, maybe iterators
impl<T> Drop for SqVec<T> {
    fn drop(&mut self) {
        if self.alloc_seg_count == 0 {
            return;
        }

        let mut seg_len: u32 = 2;
        let mut segs_left: u32 = 2;

        // seg level pointers (NonNull<NonNull<T>>)
        let mut cur_seg = self.dope;
        let last_seg = unsafe { self.dope.add(self.alloc_seg_count as usize - 1) };

        let (lseg, loff) = Self::mapping(self.len.saturating_sub(1));
        let last_nonempty_seg = unsafe { self.dope.add(lseg as usize) };

        {
            let last_el = unsafe { last_nonempty_seg.read().add(loff as usize) };
            let mut cur_item = unsafe { last_nonempty_seg.read() };
            if self.len != 0 {
                while cur_item <= last_el {
                    unsafe { cur_item.drop_in_place() };
                    cur_item = unsafe { cur_item.add(1) };
                }
            }
        }

        while cur_seg <= last_seg {
            if cur_seg < last_nonempty_seg {
                let mut cur_item = unsafe { cur_seg.read() };
                let last_item = unsafe { cur_item.add(seg_len as usize - 1) };
                while cur_item <= last_item {
                    unsafe { cur_item.drop_in_place() };
                    cur_item = unsafe { cur_item.add(1) };
                }
            }

            unsafe {
                dealloc(
                    cur_seg.read().as_ptr() as *mut u8,
                    Layout::array::<T>(seg_len as usize).unwrap_unchecked(),
                );
            }

            segs_left -= 1;
            if segs_left == 0 {
                seg_len <<= 1;
                segs_left = (3 * seg_len) >> 2;
            }

            cur_seg = unsafe { cur_seg.add(1) };
        }

        // deallocate dope vector
        unsafe {
            dealloc(
                self.dope.as_ptr() as *mut u8,
                Layout::array::<NonNull<T>>(self.dope_cap as usize).unwrap_unchecked(),
            );
        }
    }
}

unsafe impl<T: Send> Send for SqVec<T> {}
unsafe impl<T: Sync> Sync for SqVec<T> {}

impl<T> Index<u32> for SqVec<T> {
    type Output = T;

    fn index(&self, ix: u32) -> &Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_ref() }
    }
}

impl<T> IndexMut<u32> for SqVec<T> {
    fn index_mut(&mut self, ix: u32) -> &mut Self::Output {
        assert!(ix < self.len, "Index out of bounds.");
        unsafe { self.item_ptr(ix).as_mut() }
    }
}

impl<T: Debug> Debug for SqVec<T> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<T: PartialEq> PartialEq for SqVec<T> {
    fn eq(&self, other: &Self) -> bool {
        if self.len() != other.len() {
            return false;
        }
        self.iter().eq(other.iter())
    }
}

impl<T: PartialEq> PartialEq<Vec<T>> for SqVec<T> {
    fn eq(&self, other: &Vec<T>) -> bool {
        if self.len() != other.len() {
            return false;
        }
        self.iter().eq(other.iter())
    }
}

impl<T: PartialEq> PartialEq<SqVec<T>> for Vec<T> {
    fn eq(&self, other: &SqVec<T>) -> bool {
        other == self
    }
}

impl<T: Eq> Eq for SqVec<T> {}

#[cfg(test)]
mod tests {
    use super::SqVec;
    use rand::{Rng, random, rng};

    #[test]
    fn twoiter() {
        let mut sqvec = SqVec::<u32>::new();
        sqvec.push(2);
        sqvec.push(3);

        let mut iter = sqvec.iter();
        assert_eq!(iter.next(), Some(&2));
        assert_eq!(iter.next(), Some(&3));

        for _ in 0..5 {
            assert_eq!(iter.next(), None)
        }
    }

    #[test]
    fn spearfish() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        for i in 0..13 {
            sqvec.push(Box::new(i));
        }
        let mut iter = sqvec.iter();
        while let Some(item) = iter.next() {
            println!("{}", item);
        }

        iter.next();
        iter.next();
    }

    #[test]
    fn single_box_drop() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(32));
        sqvec.pop();
        assert_eq!(sqvec.iter().next(), None);
    }

    #[test]
    fn drop_len_lt_cap() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        for i in 0..200 {
            sqvec.push(Box::new(i));
        }
        for _ in 0..170 {
            sqvec.pop();
        }
    }

    #[test]
    fn drop_empty() {
        let _sqvec = SqVec::<Box<u32>>::new();
        assert_eq!(_sqvec.iter().next(), None);
    }

    #[test]
    /// randomised test of mutating api
    fn fuzz() {
        let mut rng = rng();
        let mut vec = Vec::<Box<u32>>::new();
        let mut sqvec = SqVec::<Box<u32>>::new();
        for _ in 0u32..2048 {
            let choice: u8 = rng.random_range(0..=4);
            match choice {
                0 => {
                    if !vec.is_empty() {
                        let ix = rng.random_range(0..vec.len());
                        assert_eq!(vec[ix], sqvec[ix as u32], "Not equal at index: {ix}");
                    }
                }
                1 => {
                    assert_eq!(vec.pop(), sqvec.pop());
                }
                2..=4 => {
                    let val: Box<u32> = Box::new(random());
                    vec.push(val.clone());
                    sqvec.push(val.clone());
                }
                _ => (),
            };
            assert_eq!(vec.len(), sqvec.len() as usize);
        }
    }

    #[test]
    fn lifetime() {
        let mut sqvec = SqVec::<Box<u32>>::new();
        sqvec.push(Box::new(23));
        let a = sqvec.pop().unwrap();
        drop(a);
        sqvec.push(Box::new(324));
    }

    #[test]
    fn drop_n() {
        let mut sqvec = SqVec::<u32>::new();
        for i in 0..32 {
            sqvec.push(i);
        }
    }

    #[test]
    fn iteration() {
        let mut rng = rng();
        let mut sqvec = SqVec::<u32>::new();
        let mut vec = Vec::<u32>::new();
        for _ in 0u32..600 {
            let a = rng.random();
            sqvec.push(a);
            vec.push(a);
        }
        let mut sqiter = sqvec.iter();
        let mut viter = vec.iter();
        while let (Some(a), Some(b)) = (viter.next(), sqiter.next()) {
            assert_eq!(a, b);
        }
        assert_eq!(sqiter.next(), None);
    }

    #[test]
    fn equal() {
        let mut sqvec1: SqVec<u32> = SqVec::new();
        let mut sqvec2: SqVec<u32> = SqVec::new();
        let mut vec: Vec<u32> = Vec::new();

        for i in [0, 5, 1, 56, 2] {
            sqvec1.push(i);
            sqvec2.push(i);
            vec.push(i);
        }

        assert_eq!(sqvec1, sqvec2);
        assert_eq!(sqvec1, vec);
        assert_eq!(vec, sqvec1);

        assert_ne!(sqvec!(3, 5, 2), sqvec!(2, 3));
        assert_ne!(sqvec!(33, 1, 34, 1), vec![23, 3, 5]);
        assert_ne!(sqvec!(2, 3, 6), sqvec!(2, 3, 6, 5));
    }

    #[test]
    fn iter_mut() {
        let mut rng = rng();
        let mut sqvec = SqVec::<u32>::new();
        let mut vec = Vec::<u32>::new();
        for _ in 0u32..600 {
            let r: u32 = rng.random();
            sqvec.push(r);
            vec.push(r);
        }

        let mut sqiter_mut = sqvec.iter_mut();
        let mut viter_mut = vec.iter_mut();
        while let (Some(a), Some(b)) = (viter_mut.next(), sqiter_mut.next()) {
            let r: u32 = rng.random();
            *a /= r;
            *b /= r;
            assert_eq!(a, b);
        }
    }

    // #[test]
    // fn zspearfishdei() {
    //     let sqvec: SqVec<u32> = sqvec!(2, 4, 6);
    //     let mut siter = sqvec.iter();

    //     assert_eq!(siter.next_back(), Some(&6));
    //     assert_eq!(siter.next(), Some(&2));
    //     assert_eq!(siter.next_back(), Some(&4));
    // }

    // #[test]
    // fn double_ended_iterator() {
    //     let mut rng = rng();
    //     let mut sqvec = SqVec::<u32>::new();
    //     let mut vec = Vec::<u32>::new();

    //     for _ in 0..(1 << 16) {
    //         let rnum: u32 = rng.random();
    //         sqvec.push(rnum);
    //         vec.push(rnum);
    //     }

    //     let mut sqviter = sqvec.iter();
    //     let mut viter = vec.iter();

    //     loop {
    //         let choice: bool = rng.random();

    //         if choice {
    //             match (sqviter.next(), viter.next()) {
    //                 (Some(x), Some(y)) => assert_eq!(x, y),
    //                 (None, None) => break,
    //                 _ => panic!("one iterator finished before the other"),
    //             }
    //         } else {
    //             match (sqviter.next_back(), viter.next_back()) {
    //                 (Some(x), Some(y)) => assert_eq!(x, y),
    //                 (None, None) => break,
    //                 _ => panic!("one iterator finished before the other"),
    //             }
    //         }
    //     }
    // }
}
