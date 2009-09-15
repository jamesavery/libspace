#include <assert.h>
#include <list>
#include <vector>
#include <iostream>
#include <math.h>
#include <space/storage/interval.h>

using namespace std;

#define member_local(T) template <typename Q> typename Intervals<Q>::T Intervals<Q>::
#define member(T)       template <typename Q>  T Intervals<Q>::

member_local(const_iterator) begin() const { return intervallist_t::begin(); }
member_local(const_iterator) end()   const { return intervallist_t::end(); }
member_local(iterator) begin()            { return intervallist_t::begin(); }
member_local(iterator) end()              { return intervallist_t::end(); }

member_local(overlap_status_t) overlapping(const Interval<Q>& A, const Interval<Q>& B) const {
    if(A.b < B.a) return A_STRICTLY_LESS;
    if(B.b < A.a) return B_STRICTLY_LESS;
    // Non-empty intersection
    if(A.b <= B.b) return A_LESS_OVERLAPPING;
    else return B_LESS_OVERLAPPING;
}
  
  
member(Intervals<Q>) operator *(const Intervals& B) const {
  /* Intersection of all the intervals */
  intervallist_t R;
  const_iterator iA, iB;

  for(iA=begin(), iB = B.begin();iA!=end() && iB != B.end();){
    switch(overlapping(*iA,*iB)){
    case A_STRICTLY_LESS:	// [1,3] * [4,5]
      //	cerr << "A_STRICTLY_LESS: " << (*iA) << "*" << (*iB) << " ~> dropping "<< *iA <<endl;
      iA++;
      break;
    case B_STRICTLY_LESS:    // [4,5] * [1,3]
      //	cerr << "B_STRICTLY_LESS: " << (*iA) << "*" << (*iB) << " ~> dropping " << *iB << endl;
      iB++;
      break;
    case A_LESS_OVERLAPPING: // [2,3]*[1,10] = [2,3]; iA++, [1,10]*[2,15] = [2,10]; iA++ -- altsaa A.b <= B.b
      R.push_back((*iA)*(*iB));
      //	cerr << "A_LESS_OVERLAPPING: " << (*iA) << "*" << (*iB) << " ~> " << R.back() << endl;
      iA++;
      break;
    case B_LESS_OVERLAPPING: // [1,10]*[2,3] = [2,3]; iB++, [2,15]*[1,10] = [2,10]; iB++ -- altsaa B.b <= A.b
      // Split [a,b[ i [a,B.b[, [B.b,b[; fortsaet med iA -> [B.b,b[
      R.push_back((*iA)*(*iB));
      //	cerr << "B_LESS_OVERLAPPING: " << (*iA) << "*" << (*iB) << " ~> " << R.back() << endl;
      iB++;
      break;
    }
  }
  //    cerr << "All right, we now have " << *this << endl;
  /* Throw away empty (a>=b) intervals */
  for(iterator i=R.begin();i!=R.end();++i){
    if(i->b <= i->a)
      i = R.erase(i); 
  }
  //    cerr << "all done." << endl;
  // I don't believe a merge is necessary.
  return Intervals(R);
}

// Complement
member(Intervals<Q>) operator -(const Intervals& B) const {
  // Disjoint ~> nop
  // a <= B.a && b <= B.b ~> [a;B.a[
  // a >  B.a && b >= B.b ~> [B.b;b[
  // a <  B.a && b >  B.b ~> [a;B.a[ + [B.b;b[
  // else empty?
  intervallist_t R(*this);
  iterator iA;
  const_iterator iB;

  //  cerr << endl;
  // TODO: Make linear time.
  for(iA=R.begin();iA!=R.end();iA++){
    Q a(iA->a), b(iA->b);
    for(iB=B.begin();iB!=B.end();iB++){
      if      (b <= iB->a || a >= iB->b){ 
	//	fprintf(stderr,"Disjoint: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b); // Disjoint; do nothing
      } 
      else if (a >=  iB->a && b > iB->b){
	//	fprintf(stderr,"Left overlap: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b);
	a = iB->b; // Overlap from left
      }
      else if (a < iB->a && b <= iB->b){
	//	fprintf(stderr,"Right overlap: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b);
	b = iB->a; // Overlap from right
      }
      else if (a <  iB->a && b >  iB->b){	    // Middle overlap
	//	fprintf(stderr,"Middle overlap: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b);
	iA = R.insert(++iA,interval_t(iB->b,b));
	iA--;
	b = iB->a;
      }
      else if (a>= iB->a && b<= iB->b){ // Complete cover
	//	fprintf(stderr,"Complete cover: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b);
	//	iterator i = R.erase(iA);
	//	iA = i;
	b = a;
	iB = B.end();
      } else {
	fprintf(stderr,"Unhandled: [%d;%d[\\[%d;%d[\n",a,b,iB->a,iB->b);
      }
    }
    iA->a = a;
    iA->b = b;
  }
  /* Throw away empty (a>=b) intervals */
  for(iterator i=R.begin();i!=R.end();++i){
    if(i->b <= i->a)
      i = R.erase(i); 
  }
  return Intervals(R);
}

member(Intervals<Q>) operator +(const Intervals& B) const{
  /* Union of all the intervals */
  const_iterator iA,iB;
  intervallist_t R;

  for(iA=begin(), iB = B.begin();iA!=end() && iB != B.end();){
    if(iA->b < iB->a){		// Disjoint, A smaller
      R.push_back(*iA); ++iA; 
    } else if(iB->b < iA->a){	// Disjoint, B smaller
      R.push_back(*iB); ++iB; 
    } else {			// Overlap
      R.push_back(Interval<Q>(min(iA->a,iB->a),max(iA->b,iB->b)));
      ++iA; ++iB;
    }
  }
  //  cerr << "sum: " << Intervals<Q>(R) << endl;
  // Commit excess intervals
  while(iA!=end())  { R.push_back(*iA); iA++; }
  while(iB!=B.end()){ R.push_back(*iB); iB++; }

  /* Merge overlapping intervals */
  iterator i0 = R.begin(), i1 = R.begin();
  for(++i1;i1!=R.end();){
    while(i0 != R.end() && i1 != R.end() && i0->b >= i1->a){
      const Interval<Q> I(i0->a, i1->b);
      ++i1;
      iterator i = R.erase(i0,i1);
      i0 = R.insert(i,I);
      i1 = i0;
      ++i1;
    }
    i0 = i1;
    ++i1;
  }

  return Intervals<Q>(R);
}

