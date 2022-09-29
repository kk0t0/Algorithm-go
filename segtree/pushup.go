
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
