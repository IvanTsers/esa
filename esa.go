package esa

/*
#cgo LDFLAGS: -ldivsufsort64 -L/opt/homebrew/lib
#cgo CFLAGS: -I/opt/homebrew/include
#include <divsufsort64.h>
#include <stdlib.h>
*/
import "C"
import (
	"log"
	"reflect"
	"unsafe"
)

type Stack []*Interval
type Interval struct {
	Idx int
	Lcp int
}
type Esa struct {
	T        []byte
	Sa       []int
	Lcp      []int
	Cld      []int
	KmerLen  int
	lcpCache []Minterval

	//Debug fields
	GetIntervalCallsInBuildCache int
	GetIntervalCallsInMatchPref  int
	LcpCacheHits                 int
	MatchPrefCalls               int
	MatchPrefCachedCalls         int
}
type Minterval struct {
	I, J int
	L    int
}

func (s *Stack) Top() *Interval {
	return (*s)[len(*s)-1]
}
func (s *Stack) Pop() *Interval {
	i := (*s)[len(*s)-1]
	(*s) = (*s)[0 : len(*s)-1]
	return i
}
func (s *Stack) Push(i *Interval) {
	(*s) = append(*s, i)
}
func (e *Esa) MatchPref(p []byte) *Minterval {
	//debug
	e.MatchPrefCalls += 1
	k := 0
	m := len(p)
	var parent, child *Minterval
	parent = new(Minterval)
	parent.I = 0
	parent.J = len(e.T) - 1
	for k < m {
		child = e.GetInterval(parent, p[k])
		//debug
		e.GetIntervalCallsInMatchPref += 1
		//
		if child == nil {
			parent.L = k
			return parent
		}
		l := m
		i := child.I
		j := child.J
		if i < j {
			r := 0
			if e.Lcp[i] <= e.Lcp[j+1] {
				r = e.Cld[j]
			} else {
				r = e.Cld[i]
			}
			l = min(l, e.Lcp[r])
		}
		for w := k + 1; w < l; w++ {
			if e.T[e.Sa[i]+w] != p[w] {
				child.L = w
				return child
			}
		}
		k = l
	}
	child.L = k
	return child
}
func (e *Esa) buildLcpCache() {
	kmerLen := e.KmerLen
	cacheSize := 1 << (2 * kmerLen)
	e.lcpCache = make([]Minterval, cacheSize)
	var emptyIv = Minterval{I: -1, J: -1, L: 0}
	for i := range e.lcpCache {
		e.lcpCache[i] = emptyIv
	}
	kmer := make([]byte, e.KmerLen)
	lenT := len(e.T)
	m := e.Cld[lenT-1] // left child
	rootIv := Minterval{I: 0, J: lenT - 1, L: e.Lcp[m]}
	e.walkKmerTrie(kmer, 0, rootIv)
}
func (e *Esa) walkKmerTrie(kmer []byte,
	currentK int, iv Minterval) {

	kmerLen := e.KmerLen

	if currentK < kmerLen && iv.I == -1 && iv.J == -1 {
		e.fillKmerSubtree(kmer, currentK, iv)
		return
	}
	if currentK >= kmerLen {
		e.fillKmerSubtree(kmer, currentK, iv)
		return
	}

	for code := int8(0); code < 4; code++ {
		kmer[currentK] = code2byte(code)
		parent := iv
		ij := e.GetInterval(&parent, kmer[currentK])
		e.GetIntervalCallsInBuildCache += 1

		if ij == nil {
			e.fillKmerSubtree(kmer, currentK+1, iv)
			continue
		}
		if ij.I == ij.J {
			ij.L = currentK + 1
			e.fillKmerSubtree(kmer, currentK+1, *ij)
			continue
		}
		i, j := ij.I, ij.J
		var mid int
		if e.Lcp[i] <= e.Lcp[j+1] {
			mid = e.Cld[j]
		} else {
			mid = e.Cld[i]
		}
		ij.L = e.Lcp[mid]
		if ij.L <= currentK+1 {
			e.walkKmerTrie(kmer, currentK+1, *ij)
			continue
		}
		if ij.L >= kmerLen {
			e.fillKmerSubtree(kmer, currentK+1, iv)
			continue
		}
		e.fillKmerSubtree(kmer, currentK+1, iv)

		t := currentK + 1
		nonACGT := false
		for ; t < ij.L; t++ {
			b := e.T[e.Sa[i]+t]
			if byte2code(b) < 0 {
				nonACGT = true
				break
			}
			kmer[t] = b
		}

		if nonACGT {
			e.fillKmerSubtree(kmer, t, *ij)
		} else {
			e.walkKmerTrie(kmer, t, *ij)
		}
	}
}
func (e *Esa) fillKmerSubtree(kmer []byte,
	currentK int, iv Minterval) {
	kmerLen := e.KmerLen
	if currentK < kmerLen {
		for code := int8(0); code < 4; code++ {
			kmer[currentK] = code2byte(code)
			e.fillKmerSubtree(kmer, currentK+1, iv)
		}
		return
	} else {
		idx, ok := kmer2index(kmer)
		if !ok {
			return
		}
		e.lcpCache[idx] = iv
	}
}
func (e *Esa) MatchPrefCached(p []byte) *Minterval {
	kmerLen := e.KmerLen
	e.MatchPrefCachedCalls += 1
	m := len(p)
	if m <= kmerLen {
		return e.MatchPref(p)
	}
	idx, ok := kmer2index(p[0:kmerLen])
	if !ok {
		return e.MatchPref(p)
	}
	ij := e.lcpCache[idx]
	if ij.I < 0 || (ij.I == 0 && ij.J == 0 && ij.L == 0) {
		return e.MatchPref(p)
	}
	e.LcpCacheHits += 1
	k := ij.L
	parent := &ij
	var child *Minterval
	for k < m {
		child = e.GetInterval(parent, p[k])
		//debug
		e.GetIntervalCallsInMatchPref += 1
		//
		if child == nil {
			parent.L = k
			return parent
		}
		l := m
		i := child.I
		j := child.J
		if i < j {
			r := 0
			if e.Lcp[i] <= e.Lcp[j+1] {
				r = e.Cld[j]
			} else {
				r = e.Cld[i]
			}
			l = min(l, e.Lcp[r])
		}
		for w := k + 1; w < l; w++ {
			if e.T[e.Sa[i]+w] != p[w] {
				child.L = w
				return child
			}
		}
		k = l
	}
	child.L = k
	return child
}
func Sa(t []byte) []int {
	var sa []int
	header := (*reflect.SliceHeader)(unsafe.Pointer(&t))
	ct := (*C.sauchar_t)(unsafe.Pointer(header.Data))
	n := len(t)
	csa := (*C.saidx64_t)(C.malloc(C.size_t(n * C.sizeof_saidx64_t)))
	cn := C.saidx64_t(n)
	err := int(C.divsufsort64(ct, csa, cn))
	if err != 0 {
		log.Fatalf("divsufsort failed with code %d\n", err)
	}
	header = (*reflect.SliceHeader)((unsafe.Pointer(&sa)))
	header.Cap = n
	header.Len = n
	header.Data = uintptr(unsafe.Pointer(csa))
	return sa
}
func Lcp(t []byte, sa []int) []int {
	n := len(t)
	lcp := make([]int, n)
	isa := make([]int, n)
	for i := 0; i < n; i++ {
		isa[sa[i]] = i
	}
	lcp[0] = -1
	l := 0
	for i := 0; i < n; i++ {
		j := isa[i]
		if j == 0 {
			continue
		}
		k := sa[j-1]
		for k+l < n && i+l < n && t[k+l] == t[i+l] {
			l++
		}
		lcp[j] = l
		l -= 1
		if l < 0 {
			l = 0
		}
	}
	return lcp
}
func Cld(lcp []int) []int {
	var cld []int
	lcp = append(lcp, -1)
	n := len(lcp) - 1
	cld = make([]int, n+1)
	cld[0] = n
	stack := new(Stack)
	iv := newInterval(0, -1)
	stack.Push(iv)
	for i := 1; i <= n; i++ {
		top := stack.Top()
		for lcp[i] < top.Lcp {
			last := stack.Pop()
			top = stack.Top()
			for top.Lcp == last.Lcp {
				cld[top.Idx] = last.Idx
				last = stack.Pop()
				top = stack.Top()
			}
			top = stack.Top()
			if lcp[i] < top.Lcp {
				cld[top.Idx] = last.Idx
			} else {
				cld[i-1] = last.Idx
			}
		}
		iv = newInterval(i, lcp[i])
		stack.Push(iv)
	}
	lcp = lcp[:len(lcp)-1]
	return cld
}
func newInterval(i, l int) *Interval {
	iv := new(Interval)
	iv.Idx = i
	iv.Lcp = l
	return iv
}
func MakeEsa(t []byte) *Esa {
	esa := new(Esa)
	esa.T = t
	esa.T = append(esa.T, 0)
	esa.Sa = Sa(esa.T)
	esa.Lcp = Lcp(esa.T, esa.Sa)
	esa.Lcp = append(esa.Lcp, -1)
	esa.Cld = Cld(esa.Lcp)
	esa.KmerLen = 6
	esa.buildLcpCache()
	return esa
}
func min(i, j int) int {
	if i < j {
		return i
	}
	return j
}
func (e *Esa) GetInterval(iv *Minterval, c byte) *Minterval {
	i := iv.I
	j := iv.J
	if i == j {
		if e.T[e.Sa[i]] == c {
			return iv
		}
		return nil
	}
	m := 0
	if e.Lcp[i] <= e.Lcp[j+1] {
		m = e.Cld[j]
	} else {
		m = e.Cld[i]
	}
	l := e.Lcp[m]
	k := i
	for e.Lcp[m] == l {
		if e.T[e.Sa[k]+l] == c {
			iv.I = k
			iv.J = m - 1
			return iv
		}
		k = m
		if k == j {
			break
		}
		m = e.Cld[m]
	}
	if e.T[e.Sa[k]+l] == c {
		iv.I = k
		iv.J = j
		return iv
	}
	return nil
}

var byteCodes [256]int8

func init() {
	for i := range byteCodes {
		byteCodes[i] = -1
	}
	byteCodes['A'] = 0
	byteCodes['C'] = 1
	byteCodes['G'] = 2
	byteCodes['T'] = 3
}

func byte2code(b byte) int8 {
	return byteCodes[b]
}

var codedByte = [4]byte{'A', 'C', 'G', 'T'}

func code2byte(c int8) byte {
	return codedByte[c&3]
}
func kmer2index(kmer []byte) (int, bool) {
	kmerLen := len(kmer)
	idx := 0
	for i := 0; i < kmerLen; i++ {
		c := byte2code(kmer[i])
		if c < 0 {
			return 0, false
		}
		idx = (idx << 2) | int(c)
	}
	return idx, true
}
