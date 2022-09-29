//常用的线段树单点修改区间查询操作
//维护区间最大子段和
type T struct {
	l, r, lmax, rmax, smax, sum int
}

var tr [N << 2]T
var pos [N << 2]int

func build(p, l, r int) {
	tr[p].l, tr[p].r = l, r
	if l == r {
		pos[l] = p
		tr[p].lmax, tr[p].rmax, tr[p].smax, tr[p].sum = a[l], a[l], a[l], a[l]
		return
	}
	mid := (l + r) >> 1
	build(p<<1, l, mid)
	build(p<<1|1, mid+1, r)
	pushup(&tr[p], &tr[p<<1], &tr[p<<1|1])
}
func pushup(p, l, r *T) {
	p.sum = l.sum + r.sum
	p.lmax = max(l.lmax, l.sum+r.lmax)
	p.rmax = max(r.rmax, r.sum+l.rmax)
	p.smax = max(l.smax, r.smax)
	p.smax = max(p.smax, l.rmax+r.lmax)
}
func modify(x, d int) {
	p := pos[x]
	tr[p].lmax, tr[p].rmax, tr[p].smax, tr[p].sum = d, d, d, d
	for p >>= 1; p > 0; p >>= 1 {
		pushup(&tr[p], &tr[p<<1], &tr[p<<1|1])
	}
}
func query(p, l, r int) T {
	if tr[p].l >= l && tr[p].r <= r {
		return tr[p]
	}
	mid := (tr[p].l + tr[p].r) >> 1
	if r <= mid {
		return query(p<<1, l, r)
	} else if l > mid {
		return query(p<<1|1, l, r)
	}
	left := query(p<<1, l, r)
	right := query(p<<1|1, l, r)
	var t T
	pushup(&t, &left, &right)
	return t
}
func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}
//——————————————————————————————
//线段树与dp结合
//例题：https://codeforces.com/contest/474/problem/E
 
var mp map[int64]int
var a []int64
 
type t struct {
	l, r, v, pos int
}
 
var tr [100010 << 2]t
 
func build(p, l, r int) {
	if l == r {
		tr[p].l, tr[p].r, tr[p].v = l, l, 0
		return
	}
	tr[p].l, tr[p].r = l, r
	mid := (l + r) >> 1
	build(p<<1, l, mid)
	build(p<<1|1, mid+1, r)
	pushup(p)
}
func pushup(p int) {
	if tr[p<<1].v > tr[p<<1|1].v {
		tr[p].v = tr[p<<1].v
		tr[p].pos = tr[p<<1].pos
	} else {
		tr[p].v = tr[p<<1|1].v
		tr[p].pos = tr[p<<1|1].pos
	}
}
func modify(p, x, d, pos int) {
	if tr[p].l == tr[p].r {
		if tr[p].v < d {
			tr[p].v = d
			tr[p].pos = pos
		}
		return
	}
	mid := (tr[p].l + tr[p].r) >> 1
	if x <= mid {
		modify(p<<1, x, d, pos)
	}
	if x > mid {
		modify(p<<1|1, x, d, pos)
	}
	pushup(p)
}
func query(p, l, r int) (int, int) {
	if tr[p].l >= l && tr[p].r <= r {
		return tr[p].v, tr[p].pos
	}
	mid := (tr[p].l + tr[p].r) >> 1
	pos := -1
	mmax := 0
	if l <= mid {
		a, b := query(p<<1, l, r)
		if mmax < a {
			mmax = a
			pos = b
		}
	}
	if r > mid {
		a, b := query(p<<1|1, l, r)
		if mmax < a {
			mmax = a
			pos = b
		}
	}
	return mmax, pos
}
func find(x, d int64) (int, int) {
	l := 0
	r := len(a) - 1
	t := x - d
	for l < r {
		mid := (l + r + 1) >> 1
		if a[mid] <= t {
			l = mid
		} else {
			r = mid - 1
		}
	}
	if a[l] > t {
		l--
	}
	res1 := l + 1
	l = 0
	r = len(a) - 1
	t = x + d
	for l < r {
		mid := (l + r) >> 1
		if a[mid] >= t {
			r = mid
		} else {
			l = mid + 1
		}
	}
	if a[r] < t {
		r++
	}
	res2 := r + 1
	return res1, res2
}
func solve(_r io.Reader, _w io.Writer) {
	in := bufio.NewReader(_r)
	out := bufio.NewWriter(_w)
	defer out.Flush()
	var n int
	var d int64
	Fscan(in, &n, &d)
	h := make([]int64, n+1)
	mp = make(map[int64]int, n+1)
	var v []int64
	for i := 1; i <= n; i++ {
		Fscan(in, &h[i])
		v = append(v, h[i])
	}
	if d == 0 {
		Fprintln(out, n)
		for i := 1; i <= n; i++ {
			Fprint(out, i, " ")
		}
		return
 
	}
	sort.Slice(v, func(i, j int) bool {
		return v[i] < v[j]
	})
	f := make([]int, n+1)
	pre := make([]int, n+1)
	for i := 0; i < n; i++ {
		_, ok := mp[v[i]]
		if !ok {
			mp[v[i]] = 1
			a = append(a, v[i])
		}
	}
	m := len(a)
	build(1, 1, m)
	res := 1
	for i := 1; i <= n; i++ {
		x := h[i]
		l, r := find(x, d)
		if l >= 1 {
			v1, p1 := query(1, 1, l)
			res = max(res, v1+1)
			if f[i] < v1+1 {
				f[i] = v1 + 1
				pre[i] = p1
			}
		}
		if r <= m {
			v2, p2 := query(1, r, m)
			res = max(res, v2+1)
			if f[i] < v2+1 {
				f[i] = v2 + 1
				pre[i] = p2
			}
		}
		_, t := find(x, 0)
		modify(1, t, f[i], i)
	}
	var ans []int
	for i := 1; i <= n; i++ {
		if f[i] == res {
			x := i
			for x != -1 {
				ans = append(ans, x)
				x = pre[x]
			}
			break
		}
	}
	Fprintln(out, res)
	for i := len(ans) - 1; i >= 0; i-- {
		Fprint(out, ans[i], " ")
	}
}
func main() {
	solve(os.Stdin, os.Stdout)
}
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
//势能线段树-暴力单点修改（无意义时直接return），区间查询，
//https://www.luogu.com.cn/problem/P7492
package main

import (
	"bufio"
	//"sort"

	//"crypto/des"
	. "fmt"
	"io"

	//"io/fs"

	//"strconv"

	//"math"
	"os"
	//"sort"
)

const N int = 100010

var a [N]int64

func solve(_r io.Reader, _w io.Writer) {
	in := bufio.NewReader(_r)
	out := bufio.NewWriter(_w)
	defer out.Flush()
	var n, m int
	Fscan(in, &n, &m)
	for i := 1; i <= n; i++ {
		Fscan(in, &a[i])
	}
	build(1, 1, n)
	for i := 1; i <= m; i++ {
		var opt, l, r int
		var k int64
		Fscan(in, &opt, &l, &r)
		if opt == 1 {
			Fprintln(out, max(int64(0), query(1, l, r).smax))
		} else {
			Fscan(in, &k)
			modify(1, l, r, k)
		}
	}

}
func main() {
	solve(os.Stdin, os.Stdout)
}

type T struct {
	l, r                  int
	lmax, rmax, smax, sum int64
	all                   int64 //存储补码
}

var tr [N << 2]T
var pos [N << 2]int

func build(p, l, r int) {
	tr[p].l, tr[p].r = l, r
	if l == r {
		pos[l] = p
		tr[p].lmax, tr[p].rmax, tr[p].smax, tr[p].sum = a[l], a[l], a[l], a[l]
		tr[p].all = a[l]
		return
	}
	mid := (l + r) >> 1
	build(p<<1, l, mid)
	build(p<<1|1, mid+1, r)
	pushup(&tr[p], &tr[p<<1], &tr[p<<1|1])
}
func pushup(p, l, r *T) {
	p.all = l.all & r.all
	p.sum = l.sum + r.sum
	p.lmax = max(l.lmax, l.sum+r.lmax)
	p.rmax = max(r.rmax, r.sum+l.rmax)
	p.smax = max(l.smax, r.smax)
	p.smax = max(p.smax, l.rmax+r.lmax)
}
func modify(p, l, r int, k int64) {
	if tr[p].l >= l && tr[p].r <= r {
		if (tr[p].all | k) == tr[p].all {
			return
		}
		if tr[p].l == tr[p].r {
			tr[p].lmax |= k
			tr[p].rmax |= k
			tr[p].sum |= k
			tr[p].smax |= k
			tr[p].all |= k
			return
		}
	}
	mid := (tr[p].l + tr[p].r) >> 1
	if l <= mid {
		modify(p<<1, l, r, k)
	}
	if r > mid {
		modify(p<<1|1, l, r, k)
	}
	pushup(&tr[p], &tr[p<<1], &tr[p<<1|1])
}
func query(p, l, r int) T {
	if tr[p].l >= l && tr[p].r <= r {
		return tr[p]
	}
	mid := (tr[p].l + tr[p].r) >> 1
	if r <= mid {
		return query(p<<1, l, r)
	} else if l > mid {
		return query(p<<1|1, l, r)
	}
	left := query(p<<1, l, r)
	right := query(p<<1|1, l, r)
	var t T
	pushup(&t, &left, &right)
	return t
}
func max(a, b int64) int64 {
	if a < b {
		return b
	}
	return a
}


