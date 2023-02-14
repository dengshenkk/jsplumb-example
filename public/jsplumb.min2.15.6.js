"use strict";
(function () {
  function _computeCoefficientsForAxis(curve, axis) {
    return [
      -curve[0][axis] +
        3 * curve[1][axis] +
        -3 * curve[2][axis] +
        curve[3][axis],
      3 * curve[0][axis] - 6 * curve[1][axis] + 3 * curve[2][axis],
      -3 * curve[0][axis] + 3 * curve[1][axis],
      curve[0][axis],
    ];
  }
  function _computeCoefficients(curve) {
    return [
      _computeCoefficientsForAxis(curve, "x"),
      _computeCoefficientsForAxis(curve, "y"),
    ];
  }
  function sgn(x) {
    return x < 0 ? -1 : x > 0 ? 1 : 0;
  }
  function _cubicRoots(a, b, c, d) {
    var A = b / a;
    var B = c / a;
    var C = d / a;
    var Q = (3 * B - Math.pow(A, 2)) / 9;
    var R = (9 * A * B - 27 * C - 2 * Math.pow(A, 3)) / 54;
    var D = Math.pow(Q, 3) + Math.pow(R, 2);
    var t = [];
    if (D >= 0) {
      var S =
        sgn(R + Math.sqrt(D)) * Math.pow(Math.abs(R + Math.sqrt(D)), 1 / 3);
      var T =
        sgn(R - Math.sqrt(D)) * Math.pow(Math.abs(R - Math.sqrt(D)), 1 / 3);
      t[0] = -A / 3 + (S + T);
      t[1] = -A / 3 - (S + T) / 2;
      t[2] = -A / 3 - (S + T) / 2;
      if (Math.abs((Math.sqrt(3) * (S - T)) / 2) !== 0) {
        t[1] = -1;
        t[2] = -1;
      }
    } else {
      var th = Math.acos(R / Math.sqrt(-Math.pow(Q, 3)));
      t[0] = 2 * Math.sqrt(-Q) * Math.cos(th / 3) - A / 3;
      t[1] = 2 * Math.sqrt(-Q) * Math.cos((th + 2 * Math.PI) / 3) - A / 3;
      t[2] = 2 * Math.sqrt(-Q) * Math.cos((th + 4 * Math.PI) / 3) - A / 3;
    }
    for (var i = 0; i < 3; i++) if (t[i] < 0 || t[i] > 1) t[i] = -1;
    return t;
  }
  var root = this;
  if (typeof Math.sgn == "undefined")
    Math.sgn = function (x) {
      return x == 0 ? 0 : x > 0 ? 1 : -1;
    };
  var Vectors = {
    subtract: function (v1, v2) {
      return { x: v1.x - v2.x, y: v1.y - v2.y };
    },
    dotProduct: function (v1, v2) {
      return v1.x * v2.x + v1.y * v2.y;
    },
    square: function (v) {
      return Math.sqrt(v.x * v.x + v.y * v.y);
    },
    scale: function (v, s) {
      return { x: v.x * s, y: v.y * s };
    },
  };
  var maxRecursion = 64;
  var flatnessTolerance = Math.pow(2, -maxRecursion - 1);
  var _distanceFromCurve = function (point, curve) {
    var candidates = [];
    var w = _convertToBezier(point, curve);
    var degree = curve.length - 1;
    var higherDegree = 2 * degree - 1;
    var numSolutions = _findRoots(w, higherDegree, candidates, 0);
    var v = Vectors.subtract(point, curve[0]);
    var dist = Vectors.square(v);
    var t = 0;
    for (var i = 0; i < numSolutions; i++) {
      v = Vectors.subtract(
        point,
        _bezier(curve, degree, candidates[i], null, null)
      );
      var newDist = Vectors.square(v);
      if (newDist < dist) {
        dist = newDist;
        t = candidates[i];
      }
    }
    v = Vectors.subtract(point, curve[degree]);
    newDist = Vectors.square(v);
    if (newDist < dist) {
      dist = newDist;
      t = 1;
    }
    return { location: t, distance: dist };
  };
  var _nearestPointOnCurve = function (point, curve) {
    var td = _distanceFromCurve(point, curve);
    return {
      point: _bezier(curve, curve.length - 1, td.location, null, null),
      location: td.location,
    };
  };
  var _convertToBezier = function (point, curve) {
    var degree = curve.length - 1;
    var higherDegree = 2 * degree - 1;
    var c = [];
    var d = [];
    var cdTable = [];
    var w = [];
    var z = [
      [1, 0.6, 0.3, 0.1],
      [0.4, 0.6, 0.6, 0.4],
      [0.1, 0.3, 0.6, 1],
    ];
    for (var i = 0; i <= degree; i++) c[i] = Vectors.subtract(curve[i], point);
    for (i = 0; i <= degree - 1; i++) {
      d[i] = Vectors.subtract(curve[i + 1], curve[i]);
      d[i] = Vectors.scale(d[i], 3);
    }
    for (var row = 0; row <= degree - 1; row++)
      for (var column = 0; column <= degree; column++) {
        if (!cdTable[row]) cdTable[row] = [];
        cdTable[row][column] = Vectors.dotProduct(d[row], c[column]);
      }
    for (i = 0; i <= higherDegree; i++) {
      if (!w[i]) w[i] = [];
      w[i].y = 0;
      w[i].x = parseFloat(i) / higherDegree;
    }
    var n = degree;
    var m = degree - 1;
    for (var k = 0; k <= n + m; k++) {
      var lb = Math.max(0, k - m);
      var ub = Math.min(k, n);
      for (i = lb; i <= ub; i++) {
        var j = k - i;
        w[i + j].y += cdTable[j][i] * z[j][i];
      }
    }
    return w;
  };
  var _findRoots = function (w, degree, t, depth) {
    var left = [];
    var right = [];
    var left_t = [];
    var right_t = [];
    switch (_getCrossingCount(w, degree)) {
      case 0: {
        return 0;
      }
      case 1: {
        if (depth >= maxRecursion) {
          t[0] = (w[0].x + w[degree].x) / 2;
          return 1;
        }
        if (_isFlatEnough(w, degree)) {
          t[0] = _computeXIntercept(w, degree);
          return 1;
        }
        break;
      }
    }
    _bezier(w, degree, 0.5, left, right);
    var left_count = _findRoots(left, degree, left_t, depth + 1);
    var right_count = _findRoots(right, degree, right_t, depth + 1);
    for (var i = 0; i < left_count; i++) t[i] = left_t[i];
    for (i = 0; i < right_count; i++) t[i + left_count] = right_t[i];
    return left_count + right_count;
  };
  var _getCrossingCount = function (curve, degree) {
    var n_crossings = 0;
    var old_sign;
    var sign = (old_sign = Math.sgn(curve[0].y));
    for (var i = 1; i <= degree; i++) {
      sign = Math.sgn(curve[i].y);
      if (sign != old_sign) n_crossings++;
      old_sign = sign;
    }
    return n_crossings;
  };
  var _isFlatEnough = function (curve, degree) {
    var a = curve[0].y - curve[degree].y;
    var b = curve[degree].x - curve[0].x;
    var c = curve[0].x * curve[degree].y - curve[degree].x * curve[0].y;
    var max_distance_below;
    var max_distance_above = (max_distance_below = 0);
    for (var i = 1; i < degree; i++) {
      var value = a * curve[i].x + b * curve[i].y + c;
      if (value > max_distance_above) max_distance_above = value;
      else if (value < max_distance_below) max_distance_below = value;
    }
    var a1 = 0;
    var b1 = 1;
    var c1 = 0;
    var a2 = a;
    var b2 = b;
    var c2 = c - max_distance_above;
    var det = a1 * b2 - a2 * b1;
    var dInv = 1 / det;
    var intercept_1 = (b1 * c2 - b2 * c1) * dInv;
    a2 = a;
    b2 = b;
    c2 = c - max_distance_below;
    det = a1 * b2 - a2 * b1;
    dInv = 1 / det;
    var intercept_2 = (b1 * c2 - b2 * c1) * dInv;
    var left_intercept = Math.min(intercept_1, intercept_2);
    var right_intercept = Math.max(intercept_1, intercept_2);
    var error = right_intercept - left_intercept;
    return error < flatnessTolerance ? 1 : 0;
  };
  var _computeXIntercept = function (curve, degree) {
    var XLK = 1;
    var YLK = 0;
    var XNM = curve[degree].x - curve[0].x;
    var YNM = curve[degree].y - curve[0].y;
    var XMK = curve[0].x - 0;
    var YMK = curve[0].y - 0;
    var det = XNM * YLK - YNM * XLK;
    var detInv = 1 / det;
    var S = (XNM * YMK - YNM * XMK) * detInv;
    return 0 + XLK * S;
  };
  var _bezier = function (curve, degree, t, left, right) {
    var temp = [[]];
    for (var j = 0; j <= degree; j++) temp[0][j] = curve[j];
    for (var i = 1; i <= degree; i++)
      for (j = 0; j <= degree - i; j++) {
        if (!temp[i]) temp[i] = [];
        if (!temp[i][j]) temp[i][j] = {};
        temp[i][j].x = (1 - t) * temp[i - 1][j].x + t * temp[i - 1][j + 1].x;
        temp[i][j].y = (1 - t) * temp[i - 1][j].y + t * temp[i - 1][j + 1].y;
      }
    if (left != null) for (j = 0; j <= degree; j++) left[j] = temp[j][0];
    if (right != null)
      for (j = 0; j <= degree; j++) right[j] = temp[degree - j][j];
    return temp[degree][0];
  };
  var _curveFunctionCache = {};
  var _getCurveFunctions = function (order) {
    var fns = _curveFunctionCache[order];
    if (!fns) {
      fns = [];
      var f_term = function () {
        return function (t) {
          return Math.pow(t, order);
        };
      };
      var l_term = function () {
        return function (t) {
          return Math.pow(1 - t, order);
        };
      };
      var c_term = function (c) {
        return function (t) {
          return c;
        };
      };
      var t_term = function () {
        return function (t) {
          return t;
        };
      };
      var one_minus_t_term = function () {
        return function (t) {
          return 1 - t;
        };
      };
      var _termFunc = function (terms) {
        return function (t) {
          var p = 1;
          for (var i = 0; i < terms.length; i++) p *= terms[i](t);
          return p;
        };
      };
      fns.push(new f_term());
      for (var i$jscomp$0 = 1; i$jscomp$0 < order; i$jscomp$0++) {
        var terms = [new c_term(order)];
        for (var j = 0; j < order - i$jscomp$0; j++) terms.push(new t_term());
        for (j = 0; j < i$jscomp$0; j++) terms.push(new one_minus_t_term());
        fns.push(new _termFunc(terms));
      }
      fns.push(new l_term());
      _curveFunctionCache[order] = fns;
    }
    return fns;
  };
  var _pointOnPath = function (curve, location) {
    var cc = _getCurveFunctions(curve.length - 1);
    var _x = 0;
    var _y = 0;
    for (var i = 0; i < curve.length; i++) {
      _x += curve[i].x * cc[i](location);
      _y += curve[i].y * cc[i](location);
    }
    return { x: _x, y: _y };
  };
  var _dist = function (p1, p2) {
    return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
  };
  var _isPoint = function (curve) {
    return curve[0].x === curve[1].x && curve[0].y === curve[1].y;
  };
  var _pointAlongPath = function (curve, location, distance) {
    if (_isPoint(curve)) return { point: curve[0], location };
    var prev = _pointOnPath(curve, location);
    var tally = 0;
    var curLoc = location;
    var direction = distance > 0 ? 1 : -1;
    for (var cur = null; tally < Math.abs(distance); ) {
      curLoc += 0.005 * direction;
      cur = _pointOnPath(curve, curLoc);
      tally += _dist(cur, prev);
      prev = cur;
    }
    return { point: cur, location: curLoc };
  };
  var _length = function (curve) {
    var d = new Date().getTime();
    if (_isPoint(curve)) return 0;
    var prev = _pointOnPath(curve, 0);
    var tally = 0;
    var curLoc = 0;
    var direction = 1;
    for (var cur = null; curLoc < 1; ) {
      curLoc += 0.005 * direction;
      cur = _pointOnPath(curve, curLoc);
      tally += _dist(cur, prev);
      prev = cur;
    }
    console.log("length", new Date().getTime() - d);
    return tally;
  };
  var _pointAlongPathFrom = function (curve, location, distance) {
    return _pointAlongPath(curve, location, distance).point;
  };
  var _locationAlongPathFrom = function (curve, location, distance) {
    return _pointAlongPath(curve, location, distance).location;
  };
  var _gradientAtPoint = function (curve, location) {
    var p1 = _pointOnPath(curve, location);
    var p2 = _pointOnPath(curve.slice(0, curve.length - 1), location);
    var dy = p2.y - p1.y;
    var dx = p2.x - p1.x;
    return dy === 0 ? Infinity : Math.atan(dy / dx);
  };
  var _gradientAtPointAlongPathFrom = function (curve, location, distance) {
    var p = _pointAlongPath(curve, location, distance);
    if (p.location > 1) p.location = 1;
    if (p.location < 0) p.location = 0;
    return _gradientAtPoint(curve, p.location);
  };
  var _perpendicularToPathAt = function (curve, location, length, distance) {
    distance = distance == null ? 0 : distance;
    var p = _pointAlongPath(curve, location, distance);
    var m = _gradientAtPoint(curve, p.location);
    var _theta2 = Math.atan(-1 / m);
    var y = (length / 2) * Math.sin(_theta2);
    var x = (length / 2) * Math.cos(_theta2);
    return [
      { x: p.point.x + x, y: p.point.y + y },
      { x: p.point.x - x, y: p.point.y - y },
    ];
  };
  var _lineIntersection = function (x1, y1, x2, y2, curve) {
    var a = y2 - y1;
    var b = x1 - x2;
    var c = x1 * (y1 - y2) + y1 * (x2 - x1);
    var coeffs = _computeCoefficients(curve);
    var p = [
      a * coeffs[0][0] + b * coeffs[1][0],
      a * coeffs[0][1] + b * coeffs[1][1],
      a * coeffs[0][2] + b * coeffs[1][2],
      a * coeffs[0][3] + b * coeffs[1][3] + c,
    ];
    var r = _cubicRoots.apply(null, p);
    var intersections = [];
    if (r != null)
      for (var i = 0; i < 3; i++) {
        var t = r[i];
        var t2 = Math.pow(t, 2);
        var t3 = Math.pow(t, 3);
        var x = [
          coeffs[0][0] * t3 +
            coeffs[0][1] * t2 +
            coeffs[0][2] * t +
            coeffs[0][3],
          coeffs[1][0] * t3 +
            coeffs[1][1] * t2 +
            coeffs[1][2] * t +
            coeffs[1][3],
        ];
        if (x2 - x1 !== 0) var s = (x[0] - x1) / (x2 - x1);
        else s = (x[1] - y1) / (y2 - y1);
        if (t >= 0 && t <= 1 && s >= 0 && s <= 1) intersections.push(x);
      }
    return intersections;
  };
  var _boxIntersection = function (x, y, w, h, curve) {
    var i = [];
    i.push.apply(i, _lineIntersection(x, y, x + w, y, curve));
    i.push.apply(i, _lineIntersection(x + w, y, x + w, y + h, curve));
    i.push.apply(i, _lineIntersection(x + w, y + h, x, y + h, curve));
    i.push.apply(i, _lineIntersection(x, y + h, x, y, curve));
    return i;
  };
  var _boundingBoxIntersection = function (boundingBox, curve) {
    var i = [];
    i.push.apply(
      i,
      _lineIntersection(
        boundingBox.x,
        boundingBox.y,
        boundingBox.x + boundingBox.w,
        boundingBox.y,
        curve
      )
    );
    i.push.apply(
      i,
      _lineIntersection(
        boundingBox.x + boundingBox.w,
        boundingBox.y,
        boundingBox.x + boundingBox.w,
        boundingBox.y + boundingBox.h,
        curve
      )
    );
    i.push.apply(
      i,
      _lineIntersection(
        boundingBox.x + boundingBox.w,
        boundingBox.y + boundingBox.h,
        boundingBox.x,
        boundingBox.y + boundingBox.h,
        curve
      )
    );
    i.push.apply(
      i,
      _lineIntersection(
        boundingBox.x,
        boundingBox.y + boundingBox.h,
        boundingBox.x,
        boundingBox.y,
        curve
      )
    );
    return i;
  };
  var jsBezier = (this.jsBezier = {
    distanceFromCurve: _distanceFromCurve,
    gradientAtPoint: _gradientAtPoint,
    gradientAtPointAlongCurveFrom: _gradientAtPointAlongPathFrom,
    nearestPointOnCurve: _nearestPointOnCurve,
    pointOnCurve: _pointOnPath,
    pointAlongCurveFrom: _pointAlongPathFrom,
    perpendicularToCurveAt: _perpendicularToPathAt,
    locationAlongCurveFrom: _locationAlongPathFrom,
    getLength: _length,
    lineIntersection: _lineIntersection,
    boxIntersection: _boxIntersection,
    boundingBoxIntersection: _boundingBoxIntersection,
    version: "0.9.0",
  });
  if (typeof exports !== "undefined") exports.jsBezier = jsBezier;
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var Biltong = (root.Biltong = { version: "0.4.0" });
  if (typeof exports !== "undefined") exports.Biltong = Biltong;
  var _isa = function (a) {
    return Object.prototype.toString.call(a) === "[object Array]";
  };
  var _pointHelper = function (p1, p2, fn) {
    p1 = _isa(p1) ? p1 : [p1.x, p1.y];
    p2 = _isa(p2) ? p2 : [p2.x, p2.y];
    return fn(p1, p2);
  };
  var _gradient = (Biltong.gradient = function (p1, p2) {
    return _pointHelper(p1, p2, function (_p1, _p2) {
      if (_p2[0] == _p1[0]) return _p2[1] > _p1[1] ? Infinity : -Infinity;
      else if (_p2[1] == _p1[1]) return _p2[0] > _p1[0] ? 0 : -0;
      else return (_p2[1] - _p1[1]) / (_p2[0] - _p1[0]);
    });
  });
  var _normal = (Biltong.normal = function (p1, p2) {
    return -1 / _gradient(p1, p2);
  });
  var _lineLength = (Biltong.lineLength = function (p1, p2) {
    return _pointHelper(p1, p2, function (_p1, _p2) {
      return Math.sqrt(
        Math.pow(_p2[1] - _p1[1], 2) + Math.pow(_p2[0] - _p1[0], 2)
      );
    });
  });
  var _quadrant = (Biltong.quadrant = function (p1, p2) {
    return _pointHelper(p1, p2, function (_p1, _p2) {
      if (_p2[0] > _p1[0]) return _p2[1] > _p1[1] ? 2 : 1;
      else if (_p2[0] == _p1[0]) return _p2[1] > _p1[1] ? 2 : 1;
      else return _p2[1] > _p1[1] ? 3 : 4;
    });
  });
  var _theta = (Biltong.theta = function (p1, p2) {
    return _pointHelper(p1, p2, function (_p1, _p2) {
      var m = _gradient(_p1, _p2);
      var t = Math.atan(m);
      var s = _quadrant(_p1, _p2);
      if (s == 4 || s == 3) t += Math.PI;
      if (t < 0) t += 2 * Math.PI;
      return t;
    });
  });
  var _intersects = (Biltong.intersects = function (r1, r2) {
    var x1 = r1.x;
    var x2 = r1.x + r1.w;
    var y1 = r1.y;
    var y2 = r1.y + r1.h;
    var a1 = r2.x;
    var a2 = r2.x + r2.w;
    var b1 = r2.y;
    var b2 = r2.y + r2.h;
    return (
      (x1 <= a1 && a1 <= x2 && y1 <= b1 && b1 <= y2) ||
      (x1 <= a2 && a2 <= x2 && y1 <= b1 && b1 <= y2) ||
      (x1 <= a1 && a1 <= x2 && y1 <= b2 && b2 <= y2) ||
      (x1 <= a2 && a1 <= x2 && y1 <= b2 && b2 <= y2) ||
      (a1 <= x1 && x1 <= a2 && b1 <= y1 && y1 <= b2) ||
      (a1 <= x2 && x2 <= a2 && b1 <= y1 && y1 <= b2) ||
      (a1 <= x1 && x1 <= a2 && b1 <= y2 && y2 <= b2) ||
      (a1 <= x2 && x1 <= a2 && b1 <= y2 && y2 <= b2)
    );
  });
  var _encloses = (Biltong.encloses = function (r1, r2, allowSharedEdges) {
    var x1 = r1.x;
    var x2 = r1.x + r1.w;
    var y1 = r1.y;
    var y2 = r1.y + r1.h;
    var a1 = r2.x;
    var a2 = r2.x + r2.w;
    var b1 = r2.y;
    var b2 = r2.y + r2.h;
    var c = function (v1, v2, v3, v4) {
      return allowSharedEdges ? v1 <= v2 && v3 >= v4 : v1 < v2 && v3 > v4;
    };
    return c(x1, a1, x2, a2) && c(y1, b1, y2, b2);
  });
  var _segmentMultipliers = [null, [1, -1], [1, 1], [-1, 1], [-1, -1]];
  var _inverseSegmentMultipliers = [null, [-1, -1], [-1, 1], [1, 1], [1, -1]];
  var _pointOnLine = (Biltong.pointOnLine = function (
    fromPoint,
    toPoint,
    distance
  ) {
    var m = _gradient(fromPoint, toPoint);
    var s = _quadrant(fromPoint, toPoint);
    var segmentMultiplier =
      distance > 0 ? _segmentMultipliers[s] : _inverseSegmentMultipliers[s];
    var theta = Math.atan(m);
    var y = Math.abs(distance * Math.sin(theta)) * segmentMultiplier[1];
    var x = Math.abs(distance * Math.cos(theta)) * segmentMultiplier[0];
    return { x: fromPoint.x + x, y: fromPoint.y + y };
  });
  var _perpendicularLineTo = (Biltong.perpendicularLineTo = function (
    fromPoint,
    toPoint,
    length
  ) {
    var m = _gradient(fromPoint, toPoint);
    var theta2 = Math.atan(-1 / m);
    var y = (length / 2) * Math.sin(theta2);
    var x = (length / 2) * Math.cos(theta2);
    return [
      { x: toPoint.x + x, y: toPoint.y + y },
      { x: toPoint.x - x, y: toPoint.y - y },
    ];
  });
}).call(typeof window !== "undefined" ? window : this);
(function () {
  function _touch(
    view,
    target,
    pageX,
    pageY,
    screenX,
    screenY,
    clientX,
    clientY
  ) {
    return new Touch({
      target,
      identifier: _uuid(),
      pageX,
      pageY,
      screenX,
      screenY,
      clientX: clientX || screenX,
      clientY: clientY || screenY,
    });
  }
  function _touchList() {
    var list = [];
    Array.prototype.push.apply(list, arguments);
    list.item = function (index) {
      return this[index];
    };
    return list;
  }
  function _touchAndList(
    view,
    target,
    pageX,
    pageY,
    screenX,
    screenY,
    clientX,
    clientY
  ) {
    return _touchList(_touch.apply(null, arguments));
  }
  var root = this;
  var matchesSelector = function (el, selector, ctx) {
    ctx = ctx || el.parentNode;
    var possibles = ctx.querySelectorAll(selector);
    for (var i = 0; i < possibles.length; i++)
      if (possibles[i] === el) return true;
    return false;
  };
  var _gel = function (el) {
    return typeof el == "string" || el.constructor === String
      ? document.getElementById(el)
      : el;
  };
  var _t = function (e) {
    return e.srcElement || e.target;
  };
  var _pi = function (e, target, obj, doCompute) {
    if (!doCompute) return { path: [target], end: 1 };
    else if (typeof e.path !== "undefined" && e.path.indexOf)
      return { path: e.path, end: e.path.indexOf(obj) };
    else {
      var out = { path: [], end: -1 };
      var _one = function (el) {
        out.path.push(el);
        if (el === obj) out.end = out.path.length - 1;
        else if (el.parentNode != null) _one(el.parentNode);
      };
      _one(target);
      return out;
    }
  };
  var _d = function (l, fn) {
    var i = 0;
    for (var j = l.length; i < j; i++) if (l[i] == fn) break;
    if (i < l.length) l.splice(i, 1);
  };
  var guid = 1;
  var _store = function (obj, event, fn) {
    var g = guid++;
    obj.__ta = obj.__ta || {};
    obj.__ta[event] = obj.__ta[event] || {};
    obj.__ta[event][g] = fn;
    fn.__tauid = g;
    return g;
  };
  var _unstore = function (obj, event, fn) {
    obj.__ta && obj.__ta[event] && delete obj.__ta[event][fn.__tauid];
    if (fn.__taExtra) {
      for (var i = 0; i < fn.__taExtra.length; i++)
        _unbind(obj, fn.__taExtra[i][0], fn.__taExtra[i][1]);
      fn.__taExtra.length = 0;
    }
    fn.__taUnstore && fn.__taUnstore();
  };
  var _curryChildFilter = function (children, obj, fn, evt) {
    if (children == null) return fn;
    else {
      var c = children.split(",");
      var _fn = function (e) {
        _fn.__tauid = fn.__tauid;
        var t = _t(e);
        var target = t;
        var pathInfo = _pi(e, t, obj, children != null);
        if (pathInfo.end != -1)
          for (var p = 0; p < pathInfo.end; p++) {
            target = pathInfo.path[p];
            for (var i = 0; i < c.length; i++)
              if (matchesSelector(target, c[i], obj))
                fn.apply(target, arguments);
          }
      };
      registerExtraFunction(fn, evt, _fn);
      return _fn;
    }
  };
  var registerExtraFunction = function (fn, evt, newFn) {
    fn.__taExtra = fn.__taExtra || [];
    fn.__taExtra.push([evt, newFn]);
  };
  var DefaultHandler = function (obj, evt, fn, children) {
    if (isTouchDevice && touchMap[evt]) {
      var tfn = _curryChildFilter(children, obj, fn, touchMap[evt]);
      _bind(obj, touchMap[evt], tfn, fn);
    }
    if (evt === "focus" && obj.getAttribute("tabindex") == null)
      obj.setAttribute("tabindex", "1");
    _bind(obj, evt, _curryChildFilter(children, obj, fn, evt), fn);
  };
  var SmartClickHandler = function (obj, evt, fn, children) {
    if (obj.__taSmartClicks == null) {
      var down = function (e) {
        obj.__tad = _pageLocation(e);
      };
      var up = function (e) {
        obj.__tau = _pageLocation(e);
      };
      var click = function (e) {
        if (
          obj.__tad &&
          obj.__tau &&
          obj.__tad[0] === obj.__tau[0] &&
          obj.__tad[1] === obj.__tau[1]
        )
          for (var i = 0; i < obj.__taSmartClicks.length; i++)
            obj.__taSmartClicks[i].apply(_t(e), [e]);
      };
      DefaultHandler(obj, "mousedown", down, children);
      DefaultHandler(obj, "mouseup", up, children);
      DefaultHandler(obj, "click", click, children);
      obj.__taSmartClicks = [];
    }
    obj.__taSmartClicks.push(fn);
    fn.__taUnstore = function () {
      _d(obj.__taSmartClicks, fn);
    };
  };
  var _tapProfiles = {
    tap: { touches: 1, taps: 1 },
    dbltap: { touches: 1, taps: 2 },
    contextmenu: { touches: 2, taps: 1 },
  };
  var TapHandler = function (clickThreshold, dblClickThreshold) {
    return function (obj, evt, fn, children) {
      if (evt == "contextmenu" && isMouseDevice)
        DefaultHandler(obj, evt, fn, children);
      else {
        if (obj.__taTapHandler == null) {
          var tt = (obj.__taTapHandler = {
            tap: [],
            dbltap: [],
            contextmenu: [],
            down: false,
            taps: 0,
            downSelectors: [],
          });
          var down = function (e) {
            var target = _t(e);
            var pathInfo = _pi(e, target, obj, children != null);
            var finished = false;
            for (var p = 0; p < pathInfo.end; p++) {
              if (finished) return;
              target = pathInfo.path[p];
              for (var i = 0; i < tt.downSelectors.length; i++)
                if (
                  tt.downSelectors[i] == null ||
                  matchesSelector(target, tt.downSelectors[i], obj)
                ) {
                  tt.down = true;
                  setTimeout(clearSingle, clickThreshold);
                  setTimeout(clearDouble, dblClickThreshold);
                  finished = true;
                  break;
                }
            }
          };
          var up = function (e) {
            if (tt.down) {
              var target = _t(e);
              tt.taps++;
              var tc = _touchCount(e);
              for (var eventId in _tapProfiles)
                if (_tapProfiles.hasOwnProperty(eventId)) {
                  var p = _tapProfiles[eventId];
                  if (p.touches === tc && (p.taps === 1 || p.taps === tt.taps))
                    for (var i = 0; i < tt[eventId].length; i++) {
                      var pathInfo = _pi(
                        e,
                        target,
                        obj,
                        tt[eventId][i][1] != null
                      );
                      for (var pLoop = 0; pLoop < pathInfo.end; pLoop++) {
                        var currentTarget = pathInfo.path[pLoop];
                        if (
                          tt[eventId][i][1] == null ||
                          matchesSelector(currentTarget, tt[eventId][i][1], obj)
                        ) {
                          tt[eventId][i][0].apply(currentTarget, [e]);
                          break;
                        }
                      }
                    }
                }
            }
          };
          var clearSingle = function () {
            tt.down = false;
          };
          var clearDouble = function () {
            tt.taps = 0;
          };
          DefaultHandler(obj, "mousedown", down);
          DefaultHandler(obj, "mouseup", up);
        }
        obj.__taTapHandler.downSelectors.push(children);
        obj.__taTapHandler[evt].push([fn, children]);
        fn.__taUnstore = function () {
          _d(obj.__taTapHandler[evt], fn);
        };
      }
    };
  };
  var meeHelper = function (type, evt, obj, target) {
    for (var i in obj.__tamee[type])
      if (obj.__tamee[type].hasOwnProperty(i))
        obj.__tamee[type][i].apply(target, [evt]);
  };
  var MouseEnterExitHandler = function () {
    var activeElements = [];
    return function (obj, evt, fn, children) {
      if (!obj.__tamee) {
        obj.__tamee = { over: false, mouseenter: [], mouseexit: [] };
        var over = function (e) {
          var t = _t(e);
          if (
            (children == null && t == obj && !obj.__tamee.over) ||
            (matchesSelector(t, children, obj) &&
              (t.__tamee == null || !t.__tamee.over))
          ) {
            meeHelper("mouseenter", e, obj, t);
            t.__tamee = t.__tamee || {};
            t.__tamee.over = true;
            activeElements.push(t);
          }
        };
        var out = function (e) {
          var t = _t(e);
          for (var i = 0; i < activeElements.length; i++)
            if (
              t == activeElements[i] &&
              !matchesSelector(e.relatedTarget || e.toElement, "*", t)
            ) {
              t.__tamee.over = false;
              activeElements.splice(i, 1);
              meeHelper("mouseexit", e, obj, t);
            }
        };
        _bind(
          obj,
          "mouseover",
          _curryChildFilter(children, obj, over, "mouseover"),
          over
        );
        _bind(
          obj,
          "mouseout",
          _curryChildFilter(children, obj, out, "mouseout"),
          out
        );
      }
      fn.__taUnstore = function () {
        delete obj.__tamee[evt][fn.__tauid];
      };
      _store(obj, evt, fn);
      obj.__tamee[evt][fn.__tauid] = fn;
    };
  };
  var isTouchDevice =
    "ontouchstart" in document.documentElement || navigator.maxTouchPoints;
  var isMouseDevice = "onmousedown" in document.documentElement;
  var touchMap = {
    mousedown: "touchstart",
    mouseup: "touchend",
    mousemove: "touchmove",
  };
  var touchstart = "touchstart";
  var touchend = "touchend";
  var touchmove = "touchmove";
  var iev = (function () {
    var rv = -1;
    if (navigator.appName == "Microsoft Internet Explorer") {
      var ua = navigator.userAgent;
      var re = new RegExp("MSIE ([0-9]{1,}[.0-9]{0,})");
      if (re.exec(ua) != null) rv = parseFloat(RegExp.$1);
    }
    return rv;
  })();
  var isIELT9 = iev > -1 && iev < 9;
  var _genLoc = function (e, prefix) {
    if (e == null) return [0, 0];
    var ts = _touches(e);
    var t = _getTouch(ts, 0);
    return [t[prefix + "X"], t[prefix + "Y"]];
  };
  var _pageLocation = function (e) {
    if (e == null) return [0, 0];
    if (isIELT9)
      return [
        e.clientX + document.documentElement.scrollLeft,
        e.clientY + document.documentElement.scrollTop,
      ];
    else return _genLoc(e, "page");
  };
  var _screenLocation = function (e) {
    return _genLoc(e, "screen");
  };
  var _clientLocation = function (e) {
    return _genLoc(e, "client");
  };
  var _getTouch = function (touches, idx) {
    return touches.item ? touches.item(idx) : touches[idx];
  };
  var _touches = function (e) {
    return e.touches && e.touches.length > 0
      ? e.touches
      : e.changedTouches && e.changedTouches.length > 0
      ? e.changedTouches
      : e.targetTouches && e.targetTouches.length > 0
      ? e.targetTouches
      : [e];
  };
  var _touchCount = function (e) {
    return _touches(e).length;
  };
  var _bind = function (obj, type, fn, originalFn) {
    _store(obj, type, fn);
    originalFn.__tauid = fn.__tauid;
    if (obj.addEventListener) obj.addEventListener(type, fn, false);
    else if (obj.attachEvent) {
      var key = type + fn.__tauid;
      obj["e" + key] = fn;
      obj[key] = function () {
        obj["e" + key] && obj["e" + key](window.event);
      };
      obj.attachEvent("on" + type, obj[key]);
    }
  };
  var _unbind = function (obj, type, fn) {
    if (fn == null) return;
    _each(obj, function () {
      var _el = _gel(this);
      _unstore(_el, type, fn);
      if (fn.__tauid != null)
        if (_el.removeEventListener) {
          _el.removeEventListener(type, fn, false);
          if (isTouchDevice && touchMap[type])
            _el.removeEventListener(touchMap[type], fn, false);
        } else if (this.detachEvent) {
          var key = type + fn.__tauid;
          _el[key] && _el.detachEvent("on" + type, _el[key]);
          _el[key] = null;
          _el["e" + key] = null;
        }
      if (fn.__taTouchProxy)
        _unbind(obj, fn.__taTouchProxy[1], fn.__taTouchProxy[0]);
    });
  };
  var _each = function (obj, fn) {
    if (obj == null) return;
    obj =
      typeof Window !== "undefined" &&
      typeof obj.top !== "unknown" &&
      obj == obj.top
        ? [obj]
        : typeof obj !== "string" && obj.tagName == null && obj.length != null
        ? obj
        : typeof obj === "string"
        ? document.querySelectorAll(obj)
        : [obj];
    for (var i = 0; i < obj.length; i++) fn.apply(obj[i]);
  };
  var _uuid = function () {
    return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(
      /[xy]/g,
      function (c) {
        var r = (Math.random() * 16) | 0;
        var v = c == "x" ? r : (r & 3) | 8;
        return v.toString(16);
      }
    );
  };
  root.Mottle = function (params) {
    params = params || {};
    var clickThreshold = params.clickThreshold || 250;
    var dblClickThreshold = params.dblClickThreshold || 450;
    var mouseEnterExitHandler = new MouseEnterExitHandler();
    var tapHandler = new TapHandler(clickThreshold, dblClickThreshold);
    var _smartClicks = params.smartClicks;
    var _doBind = function (obj, evt, fn, children) {
      if (fn == null) return;
      _each(obj, function () {
        var _el = _gel(this);
        if (_smartClicks && evt === "click")
          SmartClickHandler(_el, evt, fn, children);
        else if (evt === "tap" || evt === "dbltap" || evt === "contextmenu")
          tapHandler(_el, evt, fn, children);
        else if (evt === "mouseenter" || evt == "mouseexit")
          mouseEnterExitHandler(_el, evt, fn, children);
        else DefaultHandler(_el, evt, fn, children);
      });
    };
    this.remove = function (el) {
      _each(el, function () {
        var _el = _gel(this);
        if (_el.__ta)
          for (var evt in _el.__ta)
            if (_el.__ta.hasOwnProperty(evt))
              for (var h in _el.__ta[evt])
                if (_el.__ta[evt].hasOwnProperty(h))
                  _unbind(_el, evt, _el.__ta[evt][h]);
        _el.parentNode && _el.parentNode.removeChild(_el);
      });
      return this;
    };
    this.on = function (el, event, children, fn) {
      var _el = arguments[0];
      var _c = arguments.length == 4 ? arguments[2] : null;
      var _e = arguments[1];
      var _f = arguments[arguments.length - 1];
      _doBind(_el, _e, _f, _c);
      return this;
    };
    this.off = function (el, event, fn) {
      _unbind(el, event, fn);
      return this;
    };
    this.trigger = function (el, event, originalEvent, payload) {
      var originalIsMouse =
        isMouseDevice &&
        (typeof MouseEvent === "undefined" ||
          originalEvent == null ||
          originalEvent.constructor === MouseEvent);
      var eventToBind =
        isTouchDevice && !isMouseDevice && touchMap[event]
          ? touchMap[event]
          : event;
      var bindingAMouseEvent = !(
        isTouchDevice &&
        !isMouseDevice &&
        touchMap[event]
      );
      var pl = _pageLocation(originalEvent);
      var sl = _screenLocation(originalEvent);
      var cl = _clientLocation(originalEvent);
      _each(el, function () {
        var _el = _gel(this);
        originalEvent = originalEvent || {
          screenX: sl[0],
          screenY: sl[1],
          clientX: cl[0],
          clientY: cl[1],
        };
        var _decorate = function (_evt) {
          if (payload) _evt.payload = payload;
        };
        var eventGenerators = {
          TouchEvent: function (evt) {
            var touchList = _touchAndList(
              window,
              _el,
              0,
              pl[0],
              pl[1],
              sl[0],
              sl[1],
              cl[0],
              cl[1]
            );
            var init = evt.initTouchEvent || evt.initEvent;
            init(
              eventToBind,
              true,
              true,
              window,
              null,
              sl[0],
              sl[1],
              cl[0],
              cl[1],
              false,
              false,
              false,
              false,
              touchList,
              touchList,
              touchList,
              1,
              0
            );
          },
          MouseEvents: function (evt) {
            evt.initMouseEvent(
              eventToBind,
              true,
              true,
              window,
              0,
              sl[0],
              sl[1],
              cl[0],
              cl[1],
              false,
              false,
              false,
              false,
              1,
              _el
            );
          },
        };
        if (document.createEvent) {
          var ite =
            !bindingAMouseEvent &&
            !originalIsMouse &&
            isTouchDevice &&
            touchMap[event];
          var evtName = ite ? "TouchEvent" : "MouseEvents";
          var evt$jscomp$0 = document.createEvent(evtName);
          eventGenerators[evtName](evt$jscomp$0);
          _decorate(evt$jscomp$0);
          _el.dispatchEvent(evt$jscomp$0);
        } else if (document.createEventObject) {
          evt$jscomp$0 = document.createEventObject();
          evt$jscomp$0.eventType = evt$jscomp$0.eventName = eventToBind;
          evt$jscomp$0.screenX = sl[0];
          evt$jscomp$0.screenY = sl[1];
          evt$jscomp$0.clientX = cl[0];
          evt$jscomp$0.clientY = cl[1];
          _decorate(evt$jscomp$0);
          _el.fireEvent("on" + eventToBind, evt$jscomp$0);
        }
      });
      return this;
    };
  };
  root.Mottle.consume = function (e, doNotPreventDefault) {
    if (e.stopPropagation) e.stopPropagation();
    else e.returnValue = false;
    if (!doNotPreventDefault && e.preventDefault) e.preventDefault();
  };
  root.Mottle.pageLocation = _pageLocation;
  root.Mottle.setForceTouchEvents = function (value) {
    isTouchDevice = value;
  };
  root.Mottle.setForceMouseEvents = function (value) {
    isMouseDevice = value;
  };
  root.Mottle.version = "0.8.0";
  if (typeof exports !== "undefined") exports.Mottle = root.Mottle;
}).call(typeof window === "undefined" ? this : window);
(function () {
  var root = this;
  var _suggest = function (list, item, head) {
    if (list.indexOf(item) === -1) {
      head ? list.unshift(item) : list.push(item);
      return true;
    }
    return false;
  };
  var _vanquish = function (list, item) {
    var idx = list.indexOf(item);
    if (idx !== -1) list.splice(idx, 1);
  };
  var _difference = function (l1, l2) {
    var d = [];
    for (var i = 0; i < l1.length; i++)
      if (l2.indexOf(l1[i]) === -1) d.push(l1[i]);
    return d;
  };
  var _isString = function (f) {
    return f == null
      ? false
      : typeof f === "string" || f.constructor === String;
  };
  var getOffsetRect = function (elem) {
    var box = elem.getBoundingClientRect();
    var body = document.body;
    var docElem = document.documentElement;
    var scrollTop = window.pageYOffset || docElem.scrollTop || body.scrollTop;
    var scrollLeft =
      window.pageXOffset || docElem.scrollLeft || body.scrollLeft;
    var clientTop = docElem.clientTop || body.clientTop || 0;
    var clientLeft = docElem.clientLeft || body.clientLeft || 0;
    var top = box.top + scrollTop - clientTop;
    var left = box.left + scrollLeft - clientLeft;
    return { top: Math.round(top), left: Math.round(left) };
  };
  var matchesSelector = function (el, selector, ctx) {
    ctx = ctx || el.parentNode;
    var possibles = ctx.querySelectorAll(selector);
    for (var i = 0; i < possibles.length; i++)
      if (possibles[i] === el) return true;
    return false;
  };
  var findDelegateElement = function (parentElement, childElement, selector) {
    if (matchesSelector(childElement, selector, parentElement))
      return childElement;
    else
      for (
        var currentParent = childElement.parentNode;
        currentParent != null && currentParent !== parentElement;

      )
        if (matchesSelector(currentParent, selector, parentElement))
          return currentParent;
        else currentParent = currentParent.parentNode;
  };
  var findMatchingSelector = function (
    availableSelectors,
    parentElement,
    childElement
  ) {
    var el = null;
    var draggableId = parentElement.getAttribute("katavorio-draggable");
    var prefix =
      draggableId != null
        ? "[katavorio-draggable\x3d'" + draggableId + "'] "
        : "";
    for (var i = 0; i < availableSelectors.length; i++) {
      el = findDelegateElement(
        parentElement,
        childElement,
        prefix + availableSelectors[i].selector
      );
      if (el != null) {
        if (availableSelectors[i].filter) {
          var matches = matchesSelector(
            childElement,
            availableSelectors[i].filter,
            el
          );
          var exclude = availableSelectors[i].filterExclude === true;
          if ((exclude && !matches) || matches) return null;
        }
        return [availableSelectors[i], el];
      }
    }
    return null;
  };
  var iev = (function () {
    var rv = -1;
    if (navigator.appName === "Microsoft Internet Explorer") {
      var ua = navigator.userAgent;
      var re = new RegExp("MSIE ([0-9]{1,}[.0-9]{0,})");
      if (re.exec(ua) != null) rv = parseFloat(RegExp.$1);
    }
    return rv;
  })();
  var DEFAULT_GRID_X = 10;
  var DEFAULT_GRID_Y = 10;
  var isIELT9 = iev > -1 && iev < 9;
  var isIE9 = iev === 9;
  var _pl = function (e) {
    if (isIELT9)
      return [
        e.clientX + document.documentElement.scrollLeft,
        e.clientY + document.documentElement.scrollTop,
      ];
    else {
      var ts = _touches(e);
      var t = _getTouch(ts, 0);
      return isIE9
        ? [t.pageX || t.clientX, t.pageY || t.clientY]
        : [t.pageX, t.pageY];
    }
  };
  var _getTouch = function (touches, idx) {
    return touches.item ? touches.item(idx) : touches[idx];
  };
  var _touches = function (e) {
    return e.touches && e.touches.length > 0
      ? e.touches
      : e.changedTouches && e.changedTouches.length > 0
      ? e.changedTouches
      : e.targetTouches && e.targetTouches.length > 0
      ? e.targetTouches
      : [e];
  };
  var _classes = {
    delegatedDraggable: "katavorio-delegated-draggable",
    draggable: "katavorio-draggable",
    droppable: "katavorio-droppable",
    drag: "katavorio-drag",
    selected: "katavorio-drag-selected",
    active: "katavorio-drag-active",
    hover: "katavorio-drag-hover",
    noSelect: "katavorio-drag-no-select",
    ghostProxy: "katavorio-ghost-proxy",
    clonedDrag: "katavorio-clone-drag",
  };
  var _defaultScope = "katavorio-drag-scope";
  var _events = ["stop", "start", "drag", "drop", "over", "out", "beforeStart"];
  var _devNull = function () {};
  var _true = function () {
    return true;
  };
  var _foreach = function (l, fn, from) {
    for (var i = 0; i < l.length; i++) if (l[i] != from) fn(l[i]);
  };
  var _setDroppablesActive = function (dd, val, andHover, drag) {
    _foreach(dd, function (e) {
      e.setActive(val);
      if (val) e.updatePosition();
      if (andHover) e.setHover(drag, val);
    });
  };
  var _each = function (obj, fn) {
    if (obj == null) return;
    obj =
      !_isString(obj) && obj.tagName == null && obj.length != null
        ? obj
        : [obj];
    for (var i = 0; i < obj.length; i++) fn.apply(obj[i], [obj[i]]);
  };
  var _consume = function (e) {
    if (e.stopPropagation) {
      e.stopPropagation();
      e.preventDefault();
    } else e.returnValue = false;
  };
  var _defaultInputFilterSelector = "input,textarea,select,button,option";
  var _inputFilter = function (e, el, _katavorio) {
    var t = e.srcElement || e.target;
    return !matchesSelector(t, _katavorio.getInputFilterSelector(), el);
  };
  var Super = function (el, params, css, scope) {
    this.params = params || {};
    this.el = el;
    this.params.addClass(this.el, this._class);
    this.uuid = _uuid();
    var enabled = true;
    this.setEnabled = function (e) {
      enabled = e;
    };
    this.isEnabled = function () {
      return enabled;
    };
    this.toggleEnabled = function () {
      enabled = !enabled;
    };
    this.setScope = function (scopes) {
      this.scopes = scopes ? scopes.split(/\s+/) : [scope];
    };
    this.addScope = function (scopes) {
      var m = {};
      _each(this.scopes, function (s) {
        m[s] = true;
      });
      _each(scopes ? scopes.split(/\s+/) : [], function (s) {
        m[s] = true;
      });
      this.scopes = [];
      for (var i in m) this.scopes.push(i);
    };
    this.removeScope = function (scopes) {
      var m = {};
      _each(this.scopes, function (s) {
        m[s] = true;
      });
      _each(scopes ? scopes.split(/\s+/) : [], function (s) {
        delete m[s];
      });
      this.scopes = [];
      for (var i in m) this.scopes.push(i);
    };
    this.toggleScope = function (scopes) {
      var m = {};
      _each(this.scopes, function (s) {
        m[s] = true;
      });
      _each(scopes ? scopes.split(/\s+/) : [], function (s) {
        if (m[s]) delete m[s];
        else m[s] = true;
      });
      this.scopes = [];
      for (var i in m) this.scopes.push(i);
    };
    this.setScope(params.scope);
    this.k = params.katavorio;
    return params.katavorio;
  };
  var TRUE = function () {
    return true;
  };
  var FALSE = function () {
    return false;
  };
  var Drag = function (el$jscomp$0, params$jscomp$0, css, scope) {
    this._class = css.draggable;
    var k = Super.apply(this, arguments);
    this.rightButtonCanDrag = this.params.rightButtonCanDrag;
    var downAt = [0, 0];
    var posAtDown = null;
    var pagePosAtDown = null;
    var pageDelta = [0, 0];
    var moving = false;
    var initialScroll = [0, 0];
    var consumeStartEvent = this.params.consumeStartEvent !== false;
    var dragEl = this.el;
    var clone = this.params.clone;
    var scroll = this.params.scroll;
    var _multipleDrop = params$jscomp$0.multipleDrop !== false;
    var isConstrained = false;
    var elementToDrag = null;
    var availableSelectors = [];
    var activeSelectorParams = null;
    var ghostProxyParent = params$jscomp$0.ghostProxyParent;
    var currentParentPosition;
    var ghostParentPosition;
    var ghostDx;
    var ghostDy;
    if (params$jscomp$0.ghostProxy === true) var useGhostProxy = TRUE;
    else if (
      params$jscomp$0.ghostProxy &&
      typeof params$jscomp$0.ghostProxy === "function"
    )
      useGhostProxy = params$jscomp$0.ghostProxy;
    else
      useGhostProxy = function (container, dragEl) {
        if (activeSelectorParams && activeSelectorParams.useGhostProxy)
          return activeSelectorParams.useGhostProxy(container, dragEl);
        else return false;
      };
    if (params$jscomp$0.makeGhostProxy)
      var ghostProxy = params$jscomp$0.makeGhostProxy;
    else
      ghostProxy = function (el) {
        if (activeSelectorParams && activeSelectorParams.makeGhostProxy)
          return activeSelectorParams.makeGhostProxy(el);
        else return el.cloneNode(true);
      };
    if (params$jscomp$0.selector) {
      var draggableId = el$jscomp$0.getAttribute("katavorio-draggable");
      if (draggableId == null) {
        draggableId = "" + new Date().getTime();
        el$jscomp$0.setAttribute("katavorio-draggable", draggableId);
      }
      availableSelectors.push(params$jscomp$0);
    }
    var snapThreshold = params$jscomp$0.snapThreshold;
    var _snap = function (pos, gridX, gridY, thresholdX, thresholdY) {
      var _dx = Math.floor(pos[0] / gridX);
      var _dxl = gridX * _dx;
      var _dxt = _dxl + gridX;
      var _x =
        Math.abs(pos[0] - _dxl) <= thresholdX
          ? _dxl
          : Math.abs(_dxt - pos[0]) <= thresholdX
          ? _dxt
          : pos[0];
      var _dy = Math.floor(pos[1] / gridY);
      var _dyl = gridY * _dy;
      var _dyt = _dyl + gridY;
      var _y =
        Math.abs(pos[1] - _dyl) <= thresholdY
          ? _dyl
          : Math.abs(_dyt - pos[1]) <= thresholdY
          ? _dyt
          : pos[1];
      return [_x, _y];
    };
    this.posses = [];
    this.posseRoles = {};
    this.toGrid = function (pos) {
      if (this.params.grid == null) return pos;
      else {
        var tx = this.params.grid
          ? this.params.grid[0] / 2
          : snapThreshold
          ? snapThreshold
          : DEFAULT_GRID_X / 2;
        var ty = this.params.grid
          ? this.params.grid[1] / 2
          : snapThreshold
          ? snapThreshold
          : DEFAULT_GRID_Y / 2;
        return _snap(pos, this.params.grid[0], this.params.grid[1], tx, ty);
      }
    };
    this.snap = function (x, y) {
      if (dragEl == null) return;
      x = x || (this.params.grid ? this.params.grid[0] : DEFAULT_GRID_X);
      y = y || (this.params.grid ? this.params.grid[1] : DEFAULT_GRID_Y);
      var p = this.params.getPosition(dragEl);
      var tx = this.params.grid ? this.params.grid[0] / 2 : snapThreshold;
      var ty = this.params.grid ? this.params.grid[1] / 2 : snapThreshold;
      var snapped = _snap(p, x, y, tx, ty);
      this.params.setPosition(dragEl, snapped);
      return snapped;
    };
    this.setUseGhostProxy = function (val) {
      useGhostProxy = val ? TRUE : FALSE;
    };
    var constrain;
    var negativeFilter = function (pos) {
      return params$jscomp$0.allowNegative === false
        ? [Math.max(0, pos[0]), Math.max(0, pos[1])]
        : pos;
    };
    var _setConstrain = function (value) {
      constrain =
        typeof value === "function"
          ? value
          : value
          ? function (pos, dragEl, _constrainRect, _size) {
              return negativeFilter([
                Math.max(0, Math.min(_constrainRect.w - _size[0], pos[0])),
                Math.max(0, Math.min(_constrainRect.h - _size[1], pos[1])),
              ]);
            }.bind(this)
          : function (pos) {
              return negativeFilter(pos);
            };
    }.bind(this);
    _setConstrain(
      typeof this.params.constrain === "function"
        ? this.params.constrain
        : this.params.constrain || this.params.containment
    );
    this.setConstrain = function (value) {
      _setConstrain(value);
    };
    var _doConstrain = function (pos, dragEl, _constrainRect, _size) {
      if (
        activeSelectorParams != null &&
        activeSelectorParams.constrain &&
        typeof activeSelectorParams.constrain === "function"
      )
        return activeSelectorParams.constrain(
          pos,
          dragEl,
          _constrainRect,
          _size
        );
      else return constrain(pos, dragEl, _constrainRect, _size);
    };
    this.setRevert = function (fn) {
      revertFunction = fn;
    };
    if (this.params.revert) var revertFunction = this.params.revert;
    var _assignId = function (obj) {
      if (typeof obj === "function") {
        obj._katavorioId = _uuid();
        return obj._katavorioId;
      } else return obj;
    };
    var _filters = {};
    var _testFilter = function (e) {
      for (var key in _filters) {
        var f = _filters[key];
        var rv = f[0](e);
        if (f[1]) rv = !rv;
        if (!rv) return false;
      }
      return true;
    };
    var _setFilter = (this.setFilter = function (f, _exclude) {
      if (f) {
        var key = _assignId(f);
        _filters[key] = [
          function (e) {
            var t = e.srcElement || e.target;
            if (_isString(f)) var m = matchesSelector(t, f, el$jscomp$0);
            else if (typeof f === "function") m = f(e, el$jscomp$0);
            return m;
          },
          _exclude !== false,
        ];
      }
    });
    var _addFilter = (this.addFilter = _setFilter);
    var _removeFilter = (this.removeFilter = function (f) {
      var key = typeof f === "function" ? f._katavorioId : f;
      delete _filters[key];
    });
    this.clearAllFilters = function () {
      _filters = {};
    };
    this.canDrag = this.params.canDrag || _true;
    var constrainRect;
    var matchingDroppables = [];
    var intersectingDroppables = [];
    this.addSelector = function (params) {
      if (params.selector) availableSelectors.push(params);
    };
    this.downListener = function (e) {
      if (e.defaultPrevented) return;
      var isNotRightClick =
        this.rightButtonCanDrag || (e.which !== 3 && e.button !== 2);
      if (isNotRightClick && this.isEnabled() && this.canDrag()) {
        var _f = _testFilter(e) && _inputFilter(e, this.el, this.k);
        if (_f) {
          activeSelectorParams = null;
          elementToDrag = null;
          if (availableSelectors.length > 0) {
            var match = findMatchingSelector(
              availableSelectors,
              this.el,
              e.target || e.srcElement
            );
            if (match != null) {
              activeSelectorParams = match[0];
              elementToDrag = match[1];
            }
            if (elementToDrag == null) return;
          } else elementToDrag = this.el;
          if (clone) {
            dragEl = elementToDrag.cloneNode(true);
            this.params.addClass(dragEl, _classes.clonedDrag);
            dragEl.setAttribute("id", null);
            dragEl.style.position = "absolute";
            if (this.params.parent != null) {
              var p = this.params.getPosition(this.el);
              dragEl.style.left = p[0] + "px";
              dragEl.style.top = p[1] + "px";
              this.params.parent.appendChild(dragEl);
            } else {
              var b = getOffsetRect(elementToDrag);
              dragEl.style.left = b.left + "px";
              dragEl.style.top = b.top + "px";
              document.body.appendChild(dragEl);
            }
          } else dragEl = elementToDrag;
          consumeStartEvent && _consume(e);
          downAt = _pl(e);
          if (dragEl && dragEl.parentNode)
            initialScroll = [
              dragEl.parentNode.scrollLeft,
              dragEl.parentNode.scrollTop,
            ];
          this.params.bind(document, "mousemove", this.moveListener);
          this.params.bind(document, "mouseup", this.upListener);
          k.markSelection(this);
          k.markPosses(this);
          this.params.addClass(document.body, css.noSelect);
          _dispatch("beforeStart", {
            el: this.el,
            pos: posAtDown,
            e,
            drag: this,
          });
        } else if (this.params.consumeFilteredEvents) _consume(e);
      }
    }.bind(this);
    this.moveListener = function (e) {
      if (downAt) {
        if (!moving) {
          var _continue = _dispatch("start", {
            el: this.el,
            pos: posAtDown,
            e,
            drag: this,
          });
          if (_continue !== false) {
            if (!downAt) return;
            this.mark(true);
            moving = true;
          } else this.abort();
        }
        if (downAt) {
          intersectingDroppables.length = 0;
          var pos = _pl(e);
          var dx = pos[0] - downAt[0];
          var dy = pos[1] - downAt[1];
          var z = this.params.ignoreZoom ? 1 : k.getZoom();
          if (dragEl && dragEl.parentNode) {
            dx += dragEl.parentNode.scrollLeft - initialScroll[0];
            dy += dragEl.parentNode.scrollTop - initialScroll[1];
          }
          dx /= z;
          dy /= z;
          this.moveBy(dx, dy, e);
          k.updateSelection(dx, dy, this);
          k.updatePosses(dx, dy, this);
        }
      }
    }.bind(this);
    this.upListener = function (e) {
      if (downAt) {
        downAt = null;
        this.params.unbind(document, "mousemove", this.moveListener);
        this.params.unbind(document, "mouseup", this.upListener);
        this.params.removeClass(document.body, css.noSelect);
        this.unmark(e);
        k.unmarkSelection(this, e);
        k.unmarkPosses(this, e);
        this.stop(e);
        k.notifyPosseDragStop(this, e);
        moving = false;
        intersectingDroppables.length = 0;
        if (clone) {
          dragEl && dragEl.parentNode && dragEl.parentNode.removeChild(dragEl);
          dragEl = null;
        } else if (
          revertFunction &&
          revertFunction(dragEl, this.params.getPosition(dragEl)) === true
        ) {
          this.params.setPosition(dragEl, posAtDown);
          _dispatch("revert", dragEl);
        }
      }
    }.bind(this);
    this.getFilters = function () {
      return _filters;
    };
    this.abort = function () {
      if (downAt != null) this.upListener();
    };
    this.getDragElement = function (retrieveOriginalElement) {
      return retrieveOriginalElement
        ? elementToDrag || this.el
        : dragEl || this.el;
    };
    var listeners = {
      start: [],
      drag: [],
      stop: [],
      over: [],
      out: [],
      beforeStart: [],
      revert: [],
    };
    if (params$jscomp$0.events.start)
      listeners.start.push(params$jscomp$0.events.start);
    if (params$jscomp$0.events.beforeStart)
      listeners.beforeStart.push(params$jscomp$0.events.beforeStart);
    if (params$jscomp$0.events.stop)
      listeners.stop.push(params$jscomp$0.events.stop);
    if (params$jscomp$0.events.drag)
      listeners.drag.push(params$jscomp$0.events.drag);
    if (params$jscomp$0.events.revert)
      listeners.revert.push(params$jscomp$0.events.revert);
    this.on = function (evt, fn) {
      if (listeners[evt]) listeners[evt].push(fn);
    };
    this.off = function (evt, fn) {
      if (listeners[evt]) {
        var l = [];
        for (var i = 0; i < listeners[evt].length; i++)
          if (listeners[evt][i] !== fn) l.push(listeners[evt][i]);
        listeners[evt] = l;
      }
    };
    var _dispatch = function (evt, value) {
      var result = null;
      if (activeSelectorParams && activeSelectorParams[evt])
        result = activeSelectorParams[evt](value);
      else if (listeners[evt])
        for (var i = 0; i < listeners[evt].length; i++)
          try {
            var v = listeners[evt][i](value);
            if (v != null) result = v;
          } catch (e) {}
      return result;
    };
    this.notifyStart = function (e) {
      _dispatch("start", {
        el: this.el,
        pos: this.params.getPosition(dragEl),
        e,
        drag: this,
      });
    };
    this.stop = function (e, force) {
      if (force || moving) {
        var positions = [];
        var sel = k.getSelection();
        var dPos = this.params.getPosition(dragEl);
        if (sel.length > 0)
          for (var i = 0; i < sel.length; i++) {
            var p = this.params.getPosition(sel[i].el);
            positions.push([sel[i].el, { left: p[0], top: p[1] }, sel[i]]);
          }
        else positions.push([dragEl, { left: dPos[0], top: dPos[1] }, this]);
        _dispatch("stop", {
          el: dragEl,
          pos: ghostProxyOffsets || dPos,
          finalPos: dPos,
          e,
          drag: this,
          selection: positions,
        });
      }
    };
    this.mark = function (andNotify) {
      posAtDown = this.params.getPosition(dragEl);
      pagePosAtDown = this.params.getPosition(dragEl, true);
      pageDelta = [
        pagePosAtDown[0] - posAtDown[0],
        pagePosAtDown[1] - posAtDown[1],
      ];
      this.size = this.params.getSize(dragEl);
      matchingDroppables = k.getMatchingDroppables(this);
      _setDroppablesActive(matchingDroppables, true, false, this);
      this.params.addClass(dragEl, this.params.dragClass || css.drag);
      if (this.params.getConstrainingRectangle)
        var cs = this.params.getConstrainingRectangle(dragEl);
      else cs = this.params.getSize(dragEl.parentNode);
      constrainRect = { w: cs[0], h: cs[1] };
      ghostDx = 0;
      ghostDy = 0;
      if (andNotify) k.notifySelectionDragStart(this);
    };
    var ghostProxyOffsets;
    this.unmark = function (e, doNotCheckDroppables) {
      _setDroppablesActive(matchingDroppables, false, true, this);
      if (isConstrained && useGhostProxy(elementToDrag, dragEl)) {
        ghostProxyOffsets = [
          dragEl.offsetLeft - ghostDx,
          dragEl.offsetTop - ghostDy,
        ];
        dragEl.parentNode.removeChild(dragEl);
        dragEl = elementToDrag;
      } else ghostProxyOffsets = null;
      this.params.removeClass(dragEl, this.params.dragClass || css.drag);
      matchingDroppables.length = 0;
      isConstrained = false;
      if (!doNotCheckDroppables) {
        if (intersectingDroppables.length > 0 && ghostProxyOffsets)
          params$jscomp$0.setPosition(elementToDrag, ghostProxyOffsets);
        intersectingDroppables.sort(_rankSort);
        for (var i = 0; i < intersectingDroppables.length; i++) {
          var retVal = intersectingDroppables[i].drop(this, e);
          if (retVal === true) break;
        }
      }
    };
    this.moveBy = function (dx, dy, e) {
      intersectingDroppables.length = 0;
      var desiredLoc = this.toGrid([posAtDown[0] + dx, posAtDown[1] + dy]);
      var cPos = _doConstrain(desiredLoc, dragEl, constrainRect, this.size);
      if (useGhostProxy(this.el, dragEl))
        if (desiredLoc[0] !== cPos[0] || desiredLoc[1] !== cPos[1]) {
          if (!isConstrained) {
            var gp = ghostProxy(elementToDrag);
            params$jscomp$0.addClass(gp, _classes.ghostProxy);
            if (ghostProxyParent) {
              ghostProxyParent.appendChild(gp);
              currentParentPosition = params$jscomp$0.getPosition(
                elementToDrag.parentNode,
                true
              );
              ghostParentPosition = params$jscomp$0.getPosition(
                params$jscomp$0.ghostProxyParent,
                true
              );
              ghostDx = currentParentPosition[0] - ghostParentPosition[0];
              ghostDy = currentParentPosition[1] - ghostParentPosition[1];
            } else elementToDrag.parentNode.appendChild(gp);
            dragEl = gp;
            isConstrained = true;
          }
          cPos = desiredLoc;
        } else if (isConstrained) {
          dragEl.parentNode.removeChild(dragEl);
          dragEl = elementToDrag;
          isConstrained = false;
          currentParentPosition = null;
          ghostParentPosition = null;
          ghostDx = 0;
          ghostDy = 0;
        }
      var rect = { x: cPos[0], y: cPos[1], w: this.size[0], h: this.size[1] };
      var pageRect = {
        x: rect.x + pageDelta[0],
        y: rect.y + pageDelta[1],
        w: rect.w,
        h: rect.h,
      };
      var focusDropElement = null;
      this.params.setPosition(dragEl, [cPos[0] + ghostDx, cPos[1] + ghostDy]);
      for (var i = 0; i < matchingDroppables.length; i++) {
        var r2 = {
          x: matchingDroppables[i].pagePosition[0],
          y: matchingDroppables[i].pagePosition[1],
          w: matchingDroppables[i].size[0],
          h: matchingDroppables[i].size[1],
        };
        if (
          this.params.intersects(pageRect, r2) &&
          (_multipleDrop ||
            focusDropElement == null ||
            focusDropElement === matchingDroppables[i].el) &&
          matchingDroppables[i].canDrop(this)
        ) {
          if (!focusDropElement) focusDropElement = matchingDroppables[i].el;
          intersectingDroppables.push(matchingDroppables[i]);
          matchingDroppables[i].setHover(this, true, e);
        } else if (matchingDroppables[i].isHover())
          matchingDroppables[i].setHover(this, false, e);
      }
      _dispatch("drag", { el: this.el, pos: cPos, e, drag: this });
    };
    this.destroy = function () {
      this.params.unbind(this.el, "mousedown", this.downListener);
      this.params.unbind(document, "mousemove", this.moveListener);
      this.params.unbind(document, "mouseup", this.upListener);
      this.downListener = null;
      this.upListener = null;
      this.moveListener = null;
    };
    this.params.bind(this.el, "mousedown", this.downListener);
    if (this.params.handle) _setFilter(this.params.handle, false);
    else _setFilter(this.params.filter, this.params.filterExclude);
  };
  var Drop = function (el, params, css, scope) {
    this._class = css.droppable;
    this.params = params || {};
    this.rank = params.rank || 0;
    this._activeClass = this.params.activeClass || css.active;
    this._hoverClass = this.params.hoverClass || css.hover;
    Super.apply(this, arguments);
    var hover = false;
    this.allowLoopback = this.params.allowLoopback !== false;
    this.setActive = function (val) {
      this.params[val ? "addClass" : "removeClass"](this.el, this._activeClass);
    };
    this.updatePosition = function () {
      this.position = this.params.getPosition(this.el);
      this.pagePosition = this.params.getPosition(this.el, true);
      this.size = this.params.getSize(this.el);
    };
    this.canDrop =
      this.params.canDrop ||
      function (drag) {
        return true;
      };
    this.isHover = function () {
      return hover;
    };
    this.setHover = function (drag, val, e) {
      if (
        val ||
        this.el._katavorioDragHover == null ||
        this.el._katavorioDragHover === drag.el._katavorio
      ) {
        this.params[val ? "addClass" : "removeClass"](
          this.el,
          this._hoverClass
        );
        this.el._katavorioDragHover = val ? drag.el._katavorio : null;
        if (hover !== val)
          this.params.events[val ? "over" : "out"]({
            el: this.el,
            e,
            drag,
            drop: this,
          });
        hover = val;
      }
    };
    this.drop = function (drag, event) {
      return this.params.events["drop"]({ drag, e: event, drop: this });
    };
    this.destroy = function () {
      this._class = null;
      this._activeClass = null;
      this._hoverClass = null;
      hover = null;
    };
  };
  var _uuid = function () {
    return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(
      /[xy]/g,
      function (c) {
        var r = (Math.random() * 16) | 0;
        var v = c === "x" ? r : (r & 3) | 8;
        return v.toString(16);
      }
    );
  };
  var _rankSort = function (a, b) {
    return a.rank < b.rank ? 1 : a.rank > b.rank ? -1 : 0;
  };
  var _gel = function (el) {
    if (el == null) return null;
    el =
      typeof el === "string" || el.constructor === String
        ? document.getElementById(el)
        : el;
    if (el == null) return null;
    el._katavorio = el._katavorio || _uuid();
    return el;
  };
  root.Katavorio = function (katavorioParams) {
    var _selection = [];
    var _selectionMap = {};
    this._dragsByScope = {};
    this._dropsByScope = {};
    var _zoom = 1;
    var _reg = function (obj, map) {
      _each(obj, function (_obj) {
        for (var i = 0; i < _obj.scopes.length; i++) {
          map[_obj.scopes[i]] = map[_obj.scopes[i]] || [];
          map[_obj.scopes[i]].push(_obj);
        }
      });
    };
    var _unreg = function (obj, map) {
      var c = 0;
      _each(obj, function (_obj) {
        for (var i = 0; i < _obj.scopes.length; i++)
          if (map[_obj.scopes[i]]) {
            var idx = katavorioParams.indexOf(map[_obj.scopes[i]], _obj);
            if (idx !== -1) {
              map[_obj.scopes[i]].splice(idx, 1);
              c++;
            }
          }
      });
      return c > 0;
    };
    var _getMatchingDroppables = (this.getMatchingDroppables = function (drag) {
      var dd = [];
      var _m = {};
      for (var i = 0; i < drag.scopes.length; i++) {
        var _dd = this._dropsByScope[drag.scopes[i]];
        if (_dd)
          for (var j = 0; j < _dd.length; j++)
            if (
              _dd[j].canDrop(drag) &&
              !_m[_dd[j].uuid] &&
              (_dd[j].allowLoopback || _dd[j].el !== drag.el)
            ) {
              _m[_dd[j].uuid] = true;
              dd.push(_dd[j]);
            }
      }
      dd.sort(_rankSort);
      return dd;
    });
    var _prepareParams = function (p) {
      p = p || {};
      var _p = { events: {} };
      for (var i in katavorioParams) _p[i] = katavorioParams[i];
      for (i in p) _p[i] = p[i];
      for (i = 0; i < _events.length; i++)
        _p.events[_events[i]] = p[_events[i]] || _devNull;
      _p.katavorio = this;
      return _p;
    }.bind(this);
    var _mistletoe = function (existingDrag, params) {
      for (var i = 0; i < _events.length; i++)
        if (params[_events[i]]) existingDrag.on(_events[i], params[_events[i]]);
    }.bind(this);
    var _css = {};
    var overrideCss = katavorioParams.css || {};
    var _scope = katavorioParams.scope || _defaultScope;
    for (var i$jscomp$1 in _classes) _css[i$jscomp$1] = _classes[i$jscomp$1];
    for (i$jscomp$1 in overrideCss) _css[i$jscomp$1] = overrideCss[i$jscomp$1];
    var inputFilterSelector =
      katavorioParams.inputFilterSelector || _defaultInputFilterSelector;
    this.getInputFilterSelector = function () {
      return inputFilterSelector;
    };
    this.setInputFilterSelector = function (selector) {
      inputFilterSelector = selector;
      return this;
    };
    this.draggable = function (el, params) {
      var o = [];
      _each(
        el,
        function (_el) {
          _el = _gel(_el);
          if (_el != null)
            if (_el._katavorioDrag == null) {
              var p = _prepareParams(params);
              _el._katavorioDrag = new Drag(_el, p, _css, _scope);
              _reg(_el._katavorioDrag, this._dragsByScope);
              o.push(_el._katavorioDrag);
              katavorioParams.addClass(
                _el,
                p.selector ? _css.delegatedDraggable : _css.draggable
              );
            } else _mistletoe(_el._katavorioDrag, params);
        }.bind(this)
      );
      return o;
    };
    this.droppable = function (el, params) {
      var o = [];
      _each(
        el,
        function (_el) {
          _el = _gel(_el);
          if (_el != null) {
            var drop = new Drop(_el, _prepareParams(params), _css, _scope);
            _el._katavorioDrop = _el._katavorioDrop || [];
            _el._katavorioDrop.push(drop);
            _reg(drop, this._dropsByScope);
            o.push(drop);
            katavorioParams.addClass(_el, _css.droppable);
          }
        }.bind(this)
      );
      return o;
    };
    this.select = function (el) {
      _each(el, function () {
        var _el = _gel(this);
        if (_el && _el._katavorioDrag)
          if (!_selectionMap[_el._katavorio]) {
            _selection.push(_el._katavorioDrag);
            _selectionMap[_el._katavorio] = [_el, _selection.length - 1];
            katavorioParams.addClass(_el, _css.selected);
          }
      });
      return this;
    };
    this.deselect = function (el) {
      _each(el, function () {
        var _el = _gel(this);
        if (_el && _el._katavorio) {
          var e = _selectionMap[_el._katavorio];
          if (e) {
            var _s = [];
            for (var i = 0; i < _selection.length; i++)
              if (_selection[i].el !== _el) _s.push(_selection[i]);
            _selection = _s;
            delete _selectionMap[_el._katavorio];
            katavorioParams.removeClass(_el, _css.selected);
          }
        }
      });
      return this;
    };
    this.deselectAll = function () {
      for (var i in _selectionMap) {
        var d = _selectionMap[i];
        katavorioParams.removeClass(d[0], _css.selected);
      }
      _selection.length = 0;
      _selectionMap = {};
    };
    this.markSelection = function (drag) {
      _foreach(
        _selection,
        function (e) {
          e.mark();
        },
        drag
      );
    };
    this.markPosses = function (drag) {
      if (drag.posses)
        _each(drag.posses, function (p) {
          if (drag.posseRoles[p] && _posses[p])
            _foreach(
              _posses[p].members,
              function (d) {
                d.mark();
              },
              drag
            );
        });
    };
    this.unmarkSelection = function (drag, event) {
      _foreach(
        _selection,
        function (e) {
          e.unmark(event);
        },
        drag
      );
    };
    this.unmarkPosses = function (drag, event) {
      if (drag.posses)
        _each(drag.posses, function (p) {
          if (drag.posseRoles[p] && _posses[p])
            _foreach(
              _posses[p].members,
              function (d) {
                d.unmark(event, true);
              },
              drag
            );
        });
    };
    this.getSelection = function () {
      return _selection.slice(0);
    };
    this.updateSelection = function (dx, dy, drag) {
      _foreach(
        _selection,
        function (e) {
          e.moveBy(dx, dy);
        },
        drag
      );
    };
    var _posseAction = function (fn, drag) {
      if (drag.posses)
        _each(drag.posses, function (p) {
          if (drag.posseRoles[p] && _posses[p])
            _foreach(
              _posses[p].members,
              function (e) {
                fn(e);
              },
              drag
            );
        });
    };
    this.updatePosses = function (dx, dy, drag) {
      _posseAction(function (e) {
        e.moveBy(dx, dy);
      }, drag);
    };
    this.notifyPosseDragStop = function (drag, evt) {
      _posseAction(function (e) {
        e.stop(evt, true);
      }, drag);
    };
    this.notifySelectionDragStop = function (drag, evt) {
      _foreach(
        _selection,
        function (e) {
          e.stop(evt, true);
        },
        drag
      );
    };
    this.notifySelectionDragStart = function (drag, evt) {
      _foreach(
        _selection,
        function (e) {
          e.notifyStart(evt);
        },
        drag
      );
    };
    this.setZoom = function (z) {
      _zoom = z;
    };
    this.getZoom = function () {
      return _zoom;
    };
    var _scopeManip = function (kObj, scopes, map, fn) {
      _each(kObj, function (_kObj) {
        _unreg(_kObj, map);
        _kObj[fn](scopes);
        _reg(_kObj, map);
      });
    };
    _each(
      ["set", "add", "remove", "toggle"],
      function (v) {
        this[v + "Scope"] = function (el, scopes) {
          _scopeManip(
            el._katavorioDrag,
            scopes,
            this._dragsByScope,
            v + "Scope"
          );
          _scopeManip(
            el._katavorioDrop,
            scopes,
            this._dropsByScope,
            v + "Scope"
          );
        }.bind(this);
        this[v + "DragScope"] = function (el, scopes) {
          _scopeManip(
            el.constructor === Drag ? el : el._katavorioDrag,
            scopes,
            this._dragsByScope,
            v + "Scope"
          );
        }.bind(this);
        this[v + "DropScope"] = function (el, scopes) {
          _scopeManip(
            el.constructor === Drop ? el : el._katavorioDrop,
            scopes,
            this._dropsByScope,
            v + "Scope"
          );
        }.bind(this);
      }.bind(this)
    );
    this.snapToGrid = function (x, y) {
      for (var s in this._dragsByScope)
        _foreach(this._dragsByScope[s], function (d) {
          d.snap(x, y);
        });
    };
    this.getDragsForScope = function (s) {
      return this._dragsByScope[s];
    };
    this.getDropsForScope = function (s) {
      return this._dropsByScope[s];
    };
    var _destroy = function (el, type, map) {
      el = _gel(el);
      if (el[type]) {
        var selIdx = _selection.indexOf(el[type]);
        if (selIdx >= 0) _selection.splice(selIdx, 1);
        if (_unreg(el[type], map))
          _each(el[type], function (kObj) {
            kObj.destroy();
          });
        delete el[type];
      }
    };
    var _removeListener = function (el, type, evt, fn) {
      el = _gel(el);
      if (el[type]) el[type].off(evt, fn);
    };
    this.elementRemoved = function (el) {
      if (el["_katavorioDrag"]) this.destroyDraggable(el);
      if (el["_katavorioDrop"]) this.destroyDroppable(el);
    };
    this.destroyDraggable = function (el, evt, fn) {
      if (arguments.length === 1)
        _destroy(el, "_katavorioDrag", this._dragsByScope);
      else _removeListener(el, "_katavorioDrag", evt, fn);
    };
    this.destroyDroppable = function (el, evt, fn) {
      if (arguments.length === 1)
        _destroy(el, "_katavorioDrop", this._dropsByScope);
      else _removeListener(el, "_katavorioDrop", evt, fn);
    };
    this.reset = function () {
      this._dragsByScope = {};
      this._dropsByScope = {};
      _selection = [];
      _selectionMap = {};
      _posses = {};
    };
    var _posses = {};
    var _processOneSpec = function (el, _spec, dontAddExisting) {
      var posseId = _isString(_spec) ? _spec : _spec.id;
      var active = _isString(_spec) ? true : _spec.active !== false;
      var posse =
        _posses[posseId] ||
        (function () {
          var g = { name: posseId, members: [] };
          _posses[posseId] = g;
          return g;
        })();
      _each(el, function (_el) {
        if (_el._katavorioDrag) {
          if (
            dontAddExisting &&
            _el._katavorioDrag.posseRoles[posse.name] != null
          )
            return;
          _suggest(posse.members, _el._katavorioDrag);
          _suggest(_el._katavorioDrag.posses, posse.name);
          _el._katavorioDrag.posseRoles[posse.name] = active;
        }
      });
      return posse;
    };
    this.addToPosse = function (el, spec) {
      var posses = [];
      for (var i = 1; i < arguments.length; i++)
        posses.push(_processOneSpec(el, arguments[i]));
      return posses.length === 1 ? posses[0] : posses;
    };
    this.setPosse = function (el, spec) {
      var posses = [];
      for (var i$jscomp$0 = 1; i$jscomp$0 < arguments.length; i$jscomp$0++)
        posses.push(_processOneSpec(el, arguments[i$jscomp$0], true).name);
      _each(
        el,
        function (_el) {
          if (_el._katavorioDrag) {
            var diff = _difference(_el._katavorioDrag.posses, posses);
            var p = [];
            Array.prototype.push.apply(p, _el._katavorioDrag.posses);
            for (var i = 0; i < diff.length; i++)
              this.removeFromPosse(_el, diff[i]);
          }
        }.bind(this)
      );
      return posses.length === 1 ? posses[0] : posses;
    };
    this.removeFromPosse = function (el, posseId) {
      if (arguments.length < 2)
        throw new TypeError("No posse id provided for remove operation");
      for (var i = 1; i < arguments.length; i++) {
        posseId = arguments[i];
        _each(el, function (_el) {
          if (_el._katavorioDrag && _el._katavorioDrag.posses) {
            var d = _el._katavorioDrag;
            _each(posseId, function (p) {
              _vanquish(_posses[p].members, d);
              _vanquish(d.posses, p);
              delete d.posseRoles[p];
            });
          }
        });
      }
    };
    this.removeFromAllPosses = function (el) {
      _each(el, function (_el) {
        if (_el._katavorioDrag && _el._katavorioDrag.posses) {
          var d = _el._katavorioDrag;
          _each(d.posses, function (p) {
            _vanquish(_posses[p].members, d);
          });
          d.posses.length = 0;
          d.posseRoles = {};
        }
      });
    };
    this.setPosseState = function (el, posseId, state) {
      var posse = _posses[posseId];
      if (posse)
        _each(el, function (_el) {
          if (_el._katavorioDrag && _el._katavorioDrag.posses)
            _el._katavorioDrag.posseRoles[posse.name] = state;
        });
    };
  };
  root.Katavorio.version = "1.0.0";
  if (typeof exports !== "undefined") exports.Katavorio = root.Katavorio;
}).call(typeof window !== "undefined" ? window : this);
(function () {
  function isArray(a) {
    return Object.prototype.toString.call(a) === "[object Array]";
  }
  function isNumber(n) {
    return Object.prototype.toString.call(n) === "[object Number]";
  }
  function isString(s) {
    return typeof s === "string";
  }
  function isBoolean(s) {
    return typeof s === "boolean";
  }
  function isNull(s) {
    return s == null;
  }
  function isObject(o) {
    return o == null
      ? false
      : Object.prototype.toString.call(o) === "[object Object]";
  }
  function isDate(o) {
    return Object.prototype.toString.call(o) === "[object Date]";
  }
  function isFunction(o) {
    return Object.prototype.toString.call(o) === "[object Function]";
  }
  function isNamedFunction(o) {
    return isFunction(o) && o.name != null && o.name.length > 0;
  }
  function isEmpty(o) {
    for (var i in o) if (o.hasOwnProperty(i)) return false;
    return true;
  }
  function clone(a) {
    if (isString(a)) return "" + a;
    else if (isBoolean(a)) return !!a;
    else if (isDate(a)) return new Date(a.getTime());
    else if (isFunction(a)) return a;
    else if (isArray(a)) {
      var b = [];
      for (var i = 0; i < a.length; i++) b.push(clone(a[i]));
      return b;
    } else if (isObject(a)) {
      var c = {};
      for (var j in a) c[j] = clone(a[j]);
      return c;
    } else return a;
  }
  function merge(a, b, collations, overwrites) {
    var cMap = {};
    var i;
    var oMap = {};
    collations = collations || [];
    overwrites = overwrites || [];
    for (i = 0; i < collations.length; i++) cMap[collations[i]] = true;
    for (i = 0; i < overwrites.length; i++) oMap[overwrites[i]] = true;
    var c = clone(a);
    for (i in b)
      if (c[i] == null || oMap[i]) c[i] = b[i];
      else if (isString(b[i]) || isBoolean(b[i]))
        if (!cMap[i]) c[i] = b[i];
        else {
          var ar = [];
          ar.push.apply(ar, isArray(c[i]) ? c[i] : [c[i]]);
          ar.push.apply(ar, isBoolean(b[i]) ? b[i] : [b[i]]);
          c[i] = ar;
        }
      else if (isArray(b[i])) {
        ar = [];
        if (isArray(c[i])) ar.push.apply(ar, c[i]);
        ar.push.apply(ar, b[i]);
        c[i] = ar;
      } else if (isObject(b[i])) {
        if (!isObject(c[i])) c[i] = {};
        for (var j in b[i]) c[i][j] = b[i][j];
      }
    return c;
  }
  function replace(inObj, path, value) {
    if (inObj == null) return;
    var q = inObj;
    var t = q;
    path.replace(/([^\.])+/g, function (term, lc, pos, str) {
      var array = term.match(/([^\[0-9]+){1}(\[)([0-9+])/);
      var last = pos + term.length >= str.length;
      var _getArray = function () {
        return (
          t[array[1]] ||
          (function () {
            t[array[1]] = [];
            return t[array[1]];
          })()
        );
      };
      if (last)
        if (array) _getArray()[array[3]] = value;
        else t[term] = value;
      else if (array) {
        var a_1 = _getArray();
        t =
          a_1[array[3]] ||
          (function () {
            a_1[array[3]] = {};
            return a_1[array[3]];
          })();
      } else
        t =
          t[term] ||
          (function () {
            t[term] = {};
            return t[term];
          })();
      return "";
    });
    return inObj;
  }
  function functionChain(successValue, failValue, fns) {
    for (var i = 0; i < fns.length; i++) {
      var o = fns[i][0][fns[i][1]].apply(fns[i][0], fns[i][2]);
      if (o === failValue) return o;
    }
    return successValue;
  }
  function populate(model, values, functionPrefix, doNotExpandFunctions) {
    var getValue = function (fromString) {
      var matches = fromString.match(/(\${.*?})/g);
      if (matches != null)
        for (var i = 0; i < matches.length; i++) {
          var val =
            values[matches[i].substring(2, matches[i].length - 1)] || "";
          if (val != null) fromString = fromString.replace(matches[i], val);
        }
      return fromString;
    };
    var _one = function (d) {
      if (d != null)
        if (isString(d)) return getValue(d);
        else if (
          isFunction(d) &&
          !doNotExpandFunctions &&
          (functionPrefix == null ||
            (d.name || "").indexOf(functionPrefix) === 0)
        )
          return d(values);
        else if (isArray(d)) {
          var r = [];
          for (var i = 0; i < d.length; i++) r.push(_one(d[i]));
          return r;
        } else if (isObject(d)) {
          var s = {};
          for (var j in d) s[j] = _one(d[j]);
          return s;
        } else return d;
    };
    return _one(model);
  }
  function findWithFunction(a, f) {
    if (a) for (var i = 0; i < a.length; i++) if (f(a[i])) return i;
    return -1;
  }
  function removeWithFunction(a, f) {
    var idx = findWithFunction(a, f);
    if (idx > -1) a.splice(idx, 1);
    return idx !== -1;
  }
  function remove(l, v) {
    var idx = l.indexOf(v);
    if (idx > -1) l.splice(idx, 1);
    return idx !== -1;
  }
  function addWithFunction(list, item, hashFunction) {
    if (findWithFunction(list, hashFunction) === -1) list.push(item);
  }
  function addToList(map, key, value, insertAtStart) {
    var l = map[key];
    if (l == null) {
      l = [];
      map[key] = l;
    }
    l[insertAtStart ? "unshift" : "push"](value);
    return l;
  }
  function suggest(list, item, insertAtHead) {
    if (list.indexOf(item) === -1) {
      if (insertAtHead) list.unshift(item);
      else list.push(item);
      return true;
    }
    return false;
  }
  function extend(child, parent, _protoFn) {
    var i;
    parent = isArray(parent) ? parent : [parent];
    var _copyProtoChain = function (focus) {
      for (var proto = focus.__proto__; proto != null; )
        if (proto.prototype != null) {
          for (var j in proto.prototype)
            if (
              proto.prototype.hasOwnProperty(j) &&
              !child.prototype.hasOwnProperty(j)
            )
              child.prototype[j] = proto.prototype[j];
          proto = proto.prototype.__proto__;
        } else proto = null;
    };
    for (i = 0; i < parent.length; i++) {
      for (var j$jscomp$0 in parent[i].prototype)
        if (
          parent[i].prototype.hasOwnProperty(j$jscomp$0) &&
          !child.prototype.hasOwnProperty(j$jscomp$0)
        )
          child.prototype[j$jscomp$0] = parent[i].prototype[j$jscomp$0];
      _copyProtoChain(parent[i]);
    }
    var _makeFn = function (name, protoFn) {
      return function () {
        for (i = 0; i < parent.length; i++)
          if (parent[i].prototype[name])
            parent[i].prototype[name].apply(this, arguments);
        return protoFn.apply(this, arguments);
      };
    };
    var _oneSet = function (fns) {
      for (var k in fns) child.prototype[k] = _makeFn(k, fns[k]);
    };
    if (arguments.length > 2)
      for (i = 2; i < arguments.length; i++) _oneSet(arguments[i]);
    return child;
  }
  function uuid() {
    var d0 = (Math.random() * 4294967295) | 0;
    var d1 = (Math.random() * 4294967295) | 0;
    var d2 = (Math.random() * 4294967295) | 0;
    var d3 = (Math.random() * 4294967295) | 0;
    return (
      lut[d0 & 255] +
      lut[(d0 >> 8) & 255] +
      lut[(d0 >> 16) & 255] +
      lut[(d0 >> 24) & 255] +
      "-" +
      lut[d1 & 255] +
      lut[(d1 >> 8) & 255] +
      "-" +
      lut[((d1 >> 16) & 15) | 64] +
      lut[(d1 >> 24) & 255] +
      "-" +
      lut[(d2 & 63) | 128] +
      lut[(d2 >> 8) & 255] +
      "-" +
      lut[(d2 >> 16) & 255] +
      lut[(d2 >> 24) & 255] +
      lut[d3 & 255] +
      lut[(d3 >> 8) & 255] +
      lut[(d3 >> 16) & 255] +
      lut[(d3 >> 24) & 255]
    );
  }
  function fastTrim(s) {
    if (s == null) return null;
    var str = s.replace(/^\s\s*/, "");
    var ws = /\s/;
    for (var i = str.length; ws.test(str.charAt(--i)); );
    return str.slice(0, i + 1);
  }
  function each(obj, fn) {
    obj = obj.length == null || typeof obj === "string" ? [obj] : obj;
    for (var i = 0; i < obj.length; i++) fn(obj[i]);
  }
  function map$jscomp$0(obj, fn) {
    var o = [];
    for (var i = 0; i < obj.length; i++) o.push(fn(obj[i]));
    return o;
  }
  function mergeWithParents(type, map, parentAttribute) {
    parentAttribute = parentAttribute || "parent";
    var _def = function (id) {
      return id ? map[id] : null;
    };
    var _parent = function (def) {
      return def ? _def(def[parentAttribute]) : null;
    };
    var _one = function (parent, def) {
      if (parent == null) return def;
      else {
        var overrides = [
          "anchor",
          "anchors",
          "cssClass",
          "connector",
          "paintStyle",
          "hoverPaintStyle",
          "endpoint",
          "endpoints",
        ];
        if (def.mergeStrategy === "override")
          Array.prototype.push.apply(overrides, ["events", "overlays"]);
        var d_1 = merge(parent, def, [], overrides);
        return _one(_parent(parent), d_1);
      }
    };
    var _getDef = function (t) {
      if (t == null) return {};
      if (typeof t === "string") return _def(t);
      else if (t.length) {
        var done = false;
        var i = 0;
        for (var _dd = void 0; !done && i < t.length; ) {
          _dd = _getDef(t[i]);
          if (_dd) done = true;
          else i++;
        }
        return _dd;
      }
    };
    var d = _getDef(type);
    if (d) return _one(_parent(d), d);
    else return {};
  }
  function log() {
    var args = [];
    for (var _i = 0; _i < arguments.length; _i++) args[_i] = arguments[_i];
    if (jsPlumbUtil.logEnabled && typeof console !== "undefined")
      try {
        var msg = arguments[arguments.length - 1];
        console.log(msg);
      } catch (e) {}
  }
  function wrap(wrappedFunction, newFunction, returnOnThisValue) {
    return function () {
      var r = null;
      try {
        if (newFunction != null) r = newFunction.apply(this, arguments);
      } catch (e) {
        log("jsPlumb function failed : " + e);
      }
      if (
        wrappedFunction != null &&
        (returnOnThisValue == null || r !== returnOnThisValue)
      )
        try {
          r = wrappedFunction.apply(this, arguments);
        } catch (e) {
          log("wrapped function failed : " + e);
        }
      return r;
    };
  }
  function rotatePoint(point, center, rotation) {
    var radial = [point[0] - center[0], point[1] - center[1]];
    var cr = Math.cos((rotation / 360) * Math.PI * 2);
    var sr = Math.sin((rotation / 360) * Math.PI * 2);
    return [
      radial[0] * cr - radial[1] * sr + center[0],
      radial[1] * cr + radial[0] * sr + center[1],
      cr,
      sr,
    ];
  }
  function rotateAnchorOrientation(orientation, rotation) {
    var r = rotatePoint(orientation, [0, 0], rotation);
    return [Math.round(r[0]), Math.round(r[1])];
  }
  var root = this;
  root.jsPlumbUtil = root.jsPlumbUtil || {};
  var jsPlumbUtil = root.jsPlumbUtil;
  if (typeof exports !== "undefined") exports.jsPlumbUtil = jsPlumbUtil;
  jsPlumbUtil.isArray = isArray;
  jsPlumbUtil.isNumber = isNumber;
  jsPlumbUtil.isString = isString;
  jsPlumbUtil.isBoolean = isBoolean;
  jsPlumbUtil.isNull = isNull;
  jsPlumbUtil.isObject = isObject;
  jsPlumbUtil.isDate = isDate;
  jsPlumbUtil.isFunction = isFunction;
  jsPlumbUtil.isNamedFunction = isNamedFunction;
  jsPlumbUtil.isEmpty = isEmpty;
  jsPlumbUtil.clone = clone;
  jsPlumbUtil.merge = merge;
  jsPlumbUtil.replace = replace;
  jsPlumbUtil.functionChain = functionChain;
  jsPlumbUtil.populate = populate;
  jsPlumbUtil.findWithFunction = findWithFunction;
  jsPlumbUtil.removeWithFunction = removeWithFunction;
  jsPlumbUtil.remove = remove;
  jsPlumbUtil.addWithFunction = addWithFunction;
  jsPlumbUtil.addToList = addToList;
  jsPlumbUtil.suggest = suggest;
  jsPlumbUtil.extend = extend;
  var lut = [];
  for (var i$jscomp$0 = 0; i$jscomp$0 < 256; i$jscomp$0++)
    lut[i$jscomp$0] = (i$jscomp$0 < 16 ? "0" : "") + i$jscomp$0.toString(16);
  jsPlumbUtil.uuid = uuid;
  jsPlumbUtil.fastTrim = fastTrim;
  jsPlumbUtil.each = each;
  jsPlumbUtil.map = map$jscomp$0;
  jsPlumbUtil.mergeWithParents = mergeWithParents;
  jsPlumbUtil.logEnabled = true;
  jsPlumbUtil.log = log;
  jsPlumbUtil.wrap = wrap;
  var EventGenerator = (function () {
    function EventGenerator() {
      var _this = this;
      this._listeners = {};
      this.eventsSuspended = false;
      this.tick = false;
      this.eventsToDieOn = { ready: true };
      this.queue = [];
      this.bind = function (event, listener, insertAtStart) {
        var _one = function (evt) {
          addToList(_this._listeners, evt, listener, insertAtStart);
          listener.__jsPlumb = listener.__jsPlumb || {};
          listener.__jsPlumb[uuid()] = evt;
        };
        if (typeof event === "string") _one(event);
        else if (event.length != null)
          for (var i = 0; i < event.length; i++) _one(event[i]);
        return _this;
      };
      this.fire = function (event, value, originalEvent) {
        if (!this.tick) {
          this.tick = true;
          if (!this.eventsSuspended && this._listeners[event]) {
            var l = this._listeners[event].length;
            var i = 0;
            var _gone = false;
            var ret = null;
            if (
              !this.shouldFireEvent ||
              this.shouldFireEvent(event, value, originalEvent)
            )
              for (; !_gone && i < l && ret !== false; ) {
                if (this.eventsToDieOn[event])
                  this._listeners[event][i].apply(this, [value, originalEvent]);
                else
                  try {
                    ret = this._listeners[event][i].apply(this, [
                      value,
                      originalEvent,
                    ]);
                  } catch (e) {
                    log("jsPlumb: fire failed for event " + event + " : " + e);
                  }
                i++;
                if (this._listeners == null || this._listeners[event] == null)
                  _gone = true;
              }
          }
          this.tick = false;
          this._drain();
        } else this.queue.unshift(arguments);
        return this;
      };
      this._drain = function () {
        var n = _this.queue.pop();
        if (n) _this.fire.apply(_this, n);
      };
      this.unbind = function (eventOrListener, listener) {
        if (arguments.length === 0) this._listeners = {};
        else if (arguments.length === 1)
          if (typeof eventOrListener === "string")
            delete this._listeners[eventOrListener];
          else {
            if (eventOrListener.__jsPlumb) {
              var evt = void 0;
              for (var i in eventOrListener.__jsPlumb) {
                evt = eventOrListener.__jsPlumb[i];
                remove(this._listeners[evt] || [], eventOrListener);
              }
            }
          }
        else if (arguments.length === 2)
          remove(this._listeners[eventOrListener] || [], listener);
        return this;
      };
      this.getListener = function (forEvent) {
        return _this._listeners[forEvent];
      };
      this.setSuspendEvents = function (val) {
        _this.eventsSuspended = val;
      };
      this.isSuspendEvents = function () {
        return _this.eventsSuspended;
      };
      this.silently = function (fn) {
        _this.setSuspendEvents(true);
        try {
          fn();
        } catch (e) {
          log("Cannot execute silent function " + e);
        }
        _this.setSuspendEvents(false);
      };
      this.cleanupListeners = function () {
        for (var i in _this._listeners) _this._listeners[i] = null;
      };
    }
    return EventGenerator;
  })();
  jsPlumbUtil.EventGenerator = EventGenerator;
  jsPlumbUtil.rotatePoint = rotatePoint;
  jsPlumbUtil.rotateAnchorOrientation = rotateAnchorOrientation;
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  root.jsPlumbUtil.matchesSelector = function (el, selector, ctx) {
    ctx = ctx || el.parentNode;
    var possibles = ctx.querySelectorAll(selector);
    for (var i = 0; i < possibles.length; i++)
      if (possibles[i] === el) return true;
    return false;
  };
  root.jsPlumbUtil.consume = function (e, doNotPreventDefault) {
    if (e.stopPropagation) e.stopPropagation();
    else e.returnValue = false;
    if (!doNotPreventDefault && e.preventDefault) e.preventDefault();
  };
  root.jsPlumbUtil.sizeElement = function (el, x, y, w, h) {
    if (el) {
      el.style.height = h + "px";
      el.height = h;
      el.style.width = w + "px";
      el.width = w;
      el.style.left = x + "px";
      el.style.top = y + "px";
    }
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var DEFAULT_OPTIONS = {
    deriveAnchor: function (edge, index, ep, conn) {
      return {
        top: ["TopRight", "TopLeft"],
        bottom: ["BottomRight", "BottomLeft"],
      }[edge][index];
    },
  };
  var root = this;
  var ListManager = function (jsPlumbInstance, params) {
    this.count = 0;
    this.instance = jsPlumbInstance;
    this.lists = {};
    this.options = params || {};
    this.instance.addList = function (el, options) {
      return this.listManager.addList(el, options);
    };
    this.instance.removeList = function (el) {
      this.listManager.removeList(el);
    };
    this.instance.bind(
      "manageElement",
      function (p) {
        var scrollableLists = this.instance.getSelector(
          p.el,
          "[jtk-scrollable-list]"
        );
        for (var i = 0; i < scrollableLists.length; i++)
          this.addList(scrollableLists[i]);
      }.bind(this)
    );
    this.instance.bind("unmanageElement", function (p) {
      this.removeList(p.el);
    });
    this.instance.bind(
      "connection",
      function (c, evt) {
        if (evt == null) {
          this._maybeUpdateParentList(c.source);
          this._maybeUpdateParentList(c.target);
        }
      }.bind(this)
    );
  };
  root.jsPlumbListManager = ListManager;
  ListManager.prototype = {
    addList: function (el, options) {
      var dp = this.instance.extend({}, DEFAULT_OPTIONS);
      this.instance.extend(dp, this.options);
      options = this.instance.extend(dp, options || {});
      var id = [this.instance.getInstanceIndex(), this.count++].join("_");
      this.lists[id] = new List(this.instance, el, options, id);
    },
    removeList: function (el) {
      var list = this.lists[el._jsPlumbList];
      if (list) {
        list.destroy();
        delete this.lists[el._jsPlumbList];
      }
    },
    _maybeUpdateParentList: function (el) {
      var parent = el.parentNode;
      for (
        var container = this.instance.getContainer();
        parent != null && parent !== container;

      ) {
        if (
          parent._jsPlumbList != null &&
          this.lists[parent._jsPlumbList] != null
        ) {
          parent._jsPlumbScrollHandler();
          return;
        }
        parent = parent.parentNode;
      }
    },
  };
  var List = function (instance, el$jscomp$0, options, id) {
    function deriveAnchor(edge, index, ep, conn) {
      return options.anchor
        ? options.anchor
        : options.deriveAnchor(edge, index, ep, conn);
    }
    function deriveEndpoint(edge, index, ep, conn) {
      return options.deriveEndpoint
        ? options.deriveEndpoint(edge, index, ep, conn)
        : options.endpoint
        ? options.endpoint
        : ep.type;
    }
    function _maybeUpdateDraggable(el) {
      var parent = el.parentNode;
      for (
        var container = instance.getContainer();
        parent != null && parent !== container;

      ) {
        if (instance.hasClass(parent, "jtk-managed")) {
          instance.recalculateOffsets(parent);
          return;
        }
        parent = parent.parentNode;
      }
    }
    el$jscomp$0["_jsPlumbList"] = id;
    var scrollHandler = function (e) {
      var children = instance.getSelector(el$jscomp$0, ".jtk-managed");
      var elId = instance.getId(el$jscomp$0);
      for (var i = 0; i < children.length; i++) {
        if (children[i].offsetTop < el$jscomp$0.scrollTop) {
          if (!children[i]._jsPlumbProxies) {
            children[i]._jsPlumbProxies = children[i]._jsPlumbProxies || [];
            instance.select({ source: children[i] }).each(function (c) {
              instance.proxyConnection(
                c,
                0,
                el$jscomp$0,
                elId,
                function () {
                  return deriveEndpoint("top", 0, c.endpoints[0], c);
                },
                function () {
                  return deriveAnchor("top", 0, c.endpoints[0], c);
                }
              );
              children[i]._jsPlumbProxies.push([c, 0]);
            });
            instance.select({ target: children[i] }).each(function (c) {
              instance.proxyConnection(
                c,
                1,
                el$jscomp$0,
                elId,
                function () {
                  return deriveEndpoint("top", 1, c.endpoints[1], c);
                },
                function () {
                  return deriveAnchor("top", 1, c.endpoints[1], c);
                }
              );
              children[i]._jsPlumbProxies.push([c, 1]);
            });
          }
        } else if (
          children[i].offsetTop + children[i].offsetHeight >
          el$jscomp$0.scrollTop + el$jscomp$0.offsetHeight
        ) {
          if (!children[i]._jsPlumbProxies) {
            children[i]._jsPlumbProxies = children[i]._jsPlumbProxies || [];
            instance.select({ source: children[i] }).each(function (c) {
              instance.proxyConnection(
                c,
                0,
                el$jscomp$0,
                elId,
                function () {
                  return deriveEndpoint("bottom", 0, c.endpoints[0], c);
                },
                function () {
                  return deriveAnchor("bottom", 0, c.endpoints[0], c);
                }
              );
              children[i]._jsPlumbProxies.push([c, 0]);
            });
            instance.select({ target: children[i] }).each(function (c) {
              instance.proxyConnection(
                c,
                1,
                el$jscomp$0,
                elId,
                function () {
                  return deriveEndpoint("bottom", 1, c.endpoints[1], c);
                },
                function () {
                  return deriveAnchor("bottom", 1, c.endpoints[1], c);
                }
              );
              children[i]._jsPlumbProxies.push([c, 1]);
            });
          }
        } else if (children[i]._jsPlumbProxies) {
          for (var j = 0; j < children[i]._jsPlumbProxies.length; j++)
            instance.unproxyConnection(
              children[i]._jsPlumbProxies[j][0],
              children[i]._jsPlumbProxies[j][1],
              elId
            );
          delete children[i]._jsPlumbProxies;
        }
        instance.revalidate(children[i]);
      }
      _maybeUpdateDraggable(el$jscomp$0);
    };
    instance.setAttribute(el$jscomp$0, "jtk-scrollable-list", "true");
    el$jscomp$0._jsPlumbScrollHandler = scrollHandler;
    instance.on(el$jscomp$0, "scroll", scrollHandler);
    scrollHandler();
    this.destroy = function () {
      instance.off(el$jscomp$0, "scroll", scrollHandler);
      delete el$jscomp$0._jsPlumbScrollHandler;
      var children = instance.getSelector(el$jscomp$0, ".jtk-managed");
      var elId = instance.getId(el$jscomp$0);
      for (var i = 0; i < children.length; i++)
        if (children[i]._jsPlumbProxies) {
          for (var j = 0; j < children[i]._jsPlumbProxies.length; j++)
            instance.unproxyConnection(
              children[i]._jsPlumbProxies[j][0],
              children[i]._jsPlumbProxies[j][1],
              elId
            );
          delete children[i]._jsPlumbProxies;
        }
    };
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _ju = root.jsPlumbUtil;
  var _updateHoverStyle = function (component) {
    if (component._jsPlumb.paintStyle && component._jsPlumb.hoverPaintStyle) {
      var mergedHoverStyle = {};
      jsPlumb.extend(mergedHoverStyle, component._jsPlumb.paintStyle);
      jsPlumb.extend(mergedHoverStyle, component._jsPlumb.hoverPaintStyle);
      delete component._jsPlumb.hoverPaintStyle;
      if (mergedHoverStyle.gradient && component._jsPlumb.paintStyle.fill)
        delete mergedHoverStyle.gradient;
      component._jsPlumb.hoverPaintStyle = mergedHoverStyle;
    }
  };
  var events = [
    "tap",
    "dbltap",
    "click",
    "dblclick",
    "mouseover",
    "mouseout",
    "mousemove",
    "mousedown",
    "mouseup",
    "contextmenu",
  ];
  var eventFilters = { mouseout: "mouseleave", mouseexit: "mouseleave" };
  var _updateAttachedElements = function (
    component,
    state,
    timestamp,
    sourceElement
  ) {
    var affectedElements = component.getAttachedElements();
    if (affectedElements) {
      var i = 0;
      for (var j = affectedElements.length; i < j; i++)
        if (!sourceElement || sourceElement !== affectedElements[i])
          affectedElements[i].setHover(state, true, timestamp);
    }
  };
  var _splitType = function (t) {
    return t == null ? null : t.split(" ");
  };
  var _mapType = function (map, obj, typeId) {
    for (var i in obj) map[i] = typeId;
  };
  var _each = function (fn, obj) {
    obj =
      _ju.isArray(obj) || (obj.length != null && !_ju.isString(obj))
        ? obj
        : [obj];
    for (var i = 0; i < obj.length; i++)
      try {
        fn.apply(obj[i], [obj[i]]);
      } catch (e) {
        _ju.log(".each iteration failed : " + e);
      }
  };
  var _applyTypes = function (component, params, doNotRepaint) {
    if (component.getDefaultType) {
      var td = component.getTypeDescriptor();
      var map = {};
      var defType = component.getDefaultType();
      var o = _ju.merge({}, defType);
      _mapType(map, defType, "__default");
      var i = 0;
      for (var j = component._jsPlumb.types.length; i < j; i++) {
        var tid = component._jsPlumb.types[i];
        if (tid !== "__default") {
          var _t = component._jsPlumb.instance.getType(tid, td);
          if (_t != null) {
            var overrides = [
              "anchor",
              "anchors",
              "connector",
              "paintStyle",
              "hoverPaintStyle",
              "endpoint",
              "endpoints",
              "connectorOverlays",
              "connectorStyle",
              "connectorHoverStyle",
              "endpointStyle",
              "endpointHoverStyle",
            ];
            var collations = [];
            if (_t.mergeStrategy === "override")
              Array.prototype.push.apply(overrides, [
                "events",
                "overlays",
                "cssClass",
              ]);
            else collations.push("cssClass");
            o = _ju.merge(o, _t, collations, overrides);
            _mapType(map, _t, tid);
          }
        }
      }
      if (params) o = _ju.populate(o, params, "_");
      component.applyType(o, doNotRepaint, map);
      if (!doNotRepaint) component.repaint();
    }
  };
  var jsPlumbUIComponent = (root.jsPlumbUIComponent = function (params) {
    _ju.EventGenerator.apply(this, arguments);
    var self = this;
    var a = arguments;
    var idPrefix = self.idPrefix;
    var id = idPrefix + new Date().getTime();
    this._jsPlumb = {
      instance: params._jsPlumb,
      parameters: params.parameters || {},
      paintStyle: null,
      hoverPaintStyle: null,
      paintStyleInUse: null,
      hover: false,
      beforeDetach: params.beforeDetach,
      beforeDrop: params.beforeDrop,
      overlayPlacements: [],
      hoverClass: params.hoverClass || params._jsPlumb.Defaults.HoverClass,
      types: [],
      typeCache: {},
    };
    this.cacheTypeItem = function (key, item, typeId) {
      this._jsPlumb.typeCache[typeId] = this._jsPlumb.typeCache[typeId] || {};
      this._jsPlumb.typeCache[typeId][key] = item;
    };
    this.getCachedTypeItem = function (key, typeId) {
      return this._jsPlumb.typeCache[typeId]
        ? this._jsPlumb.typeCache[typeId][key]
        : null;
    };
    this.getId = function () {
      return id;
    };
    var o$jscomp$0 = params.overlays || [];
    var oo = {};
    if (this.defaultOverlayKeys) {
      for (
        var i$jscomp$0 = 0;
        i$jscomp$0 < this.defaultOverlayKeys.length;
        i$jscomp$0++
      )
        Array.prototype.push.apply(
          o$jscomp$0,
          this._jsPlumb.instance.Defaults[
            this.defaultOverlayKeys[i$jscomp$0]
          ] || []
        );
      for (i$jscomp$0 = 0; i$jscomp$0 < o$jscomp$0.length; i$jscomp$0++) {
        var fo = jsPlumb.convertToFullOverlaySpec(o$jscomp$0[i$jscomp$0]);
        oo[fo[1].id] = fo;
      }
    }
    var _defaultType = {
      overlays: oo,
      parameters: params.parameters || {},
      scope: params.scope || this._jsPlumb.instance.getDefaultScope(),
    };
    this.getDefaultType = function () {
      return _defaultType;
    };
    this.appendToDefaultType = function (obj) {
      for (var i in obj) _defaultType[i] = obj[i];
    };
    if (params.events)
      for (var evtName in params.events)
        self.bind(evtName, params.events[evtName]);
    this.clone = function () {
      var o = Object.create(this.constructor.prototype);
      this.constructor.apply(o, a);
      return o;
    }.bind(this);
    this.isDetachAllowed = function (connection) {
      var r = true;
      if (this._jsPlumb.beforeDetach)
        try {
          r = this._jsPlumb.beforeDetach(connection);
        } catch (e) {
          _ju.log("jsPlumb: beforeDetach callback failed", e);
        }
      return r;
    };
    this.isDropAllowed = function (
      sourceId,
      targetId,
      scope,
      connection,
      dropEndpoint,
      source,
      target
    ) {
      var r = this._jsPlumb.instance.checkCondition("beforeDrop", {
        sourceId,
        targetId,
        scope,
        connection,
        dropEndpoint,
        source,
        target,
      });
      if (this._jsPlumb.beforeDrop)
        try {
          r = this._jsPlumb.beforeDrop({
            sourceId,
            targetId,
            scope,
            connection,
            dropEndpoint,
            source,
            target,
          });
        } catch (e) {
          _ju.log("jsPlumb: beforeDrop callback failed", e);
        }
      return r;
    };
    var domListeners = [];
    this.setListenerComponent = function (c) {
      for (var i = 0; i < domListeners.length; i++) domListeners[i][3] = c;
    };
  });
  var _removeTypeCssHelper = function (component, typeIndex) {
    var typeId = component._jsPlumb.types[typeIndex];
    var type = component._jsPlumb.instance.getType(
      typeId,
      component.getTypeDescriptor()
    );
    if (type != null && type.cssClass && component.canvas)
      component._jsPlumb.instance.removeClass(component.canvas, type.cssClass);
  };
  _ju.extend(root.jsPlumbUIComponent, _ju.EventGenerator, {
    getParameter: function (name) {
      return this._jsPlumb.parameters[name];
    },
    setParameter: function (name, value) {
      this._jsPlumb.parameters[name] = value;
    },
    getParameters: function () {
      return this._jsPlumb.parameters;
    },
    setParameters: function (p) {
      this._jsPlumb.parameters = p;
    },
    getClass: function () {
      return jsPlumb.getClass(this.canvas);
    },
    hasClass: function (clazz) {
      return jsPlumb.hasClass(this.canvas, clazz);
    },
    addClass: function (clazz) {
      jsPlumb.addClass(this.canvas, clazz);
    },
    removeClass: function (clazz) {
      jsPlumb.removeClass(this.canvas, clazz);
    },
    updateClasses: function (classesToAdd, classesToRemove) {
      jsPlumb.updateClasses(this.canvas, classesToAdd, classesToRemove);
    },
    setType: function (typeId, params, doNotRepaint) {
      this.clearTypes();
      this._jsPlumb.types = _splitType(typeId) || [];
      _applyTypes(this, params, doNotRepaint);
    },
    getType: function () {
      return this._jsPlumb.types;
    },
    reapplyTypes: function (params, doNotRepaint) {
      _applyTypes(this, params, doNotRepaint);
    },
    hasType: function (typeId) {
      return this._jsPlumb.types.indexOf(typeId) !== -1;
    },
    addType: function (typeId, params, doNotRepaint) {
      var t = _splitType(typeId);
      var _cont = false;
      if (t != null) {
        var i = 0;
        for (var j = t.length; i < j; i++)
          if (!this.hasType(t[i])) {
            this._jsPlumb.types.push(t[i]);
            _cont = true;
          }
        if (_cont) _applyTypes(this, params, doNotRepaint);
      }
    },
    removeType: function (typeId, params, doNotRepaint) {
      var t = _splitType(typeId);
      var _cont = false;
      var _one = function (tt) {
        var idx = this._jsPlumb.types.indexOf(tt);
        if (idx !== -1) {
          _removeTypeCssHelper(this, idx);
          this._jsPlumb.types.splice(idx, 1);
          return true;
        }
        return false;
      }.bind(this);
      if (t != null) {
        var i = 0;
        for (var j = t.length; i < j; i++) _cont = _one(t[i]) || _cont;
        if (_cont) _applyTypes(this, params, doNotRepaint);
      }
    },
    clearTypes: function (params, doNotRepaint) {
      var i = this._jsPlumb.types.length;
      for (var j = 0; j < i; j++) {
        _removeTypeCssHelper(this, 0);
        this._jsPlumb.types.splice(0, 1);
      }
      _applyTypes(this, params, doNotRepaint);
    },
    toggleType: function (typeId, params, doNotRepaint) {
      var t = _splitType(typeId);
      if (t != null) {
        var i = 0;
        for (var j = t.length; i < j; i++) {
          var idx = this._jsPlumb.types.indexOf(t[i]);
          if (idx !== -1) {
            _removeTypeCssHelper(this, idx);
            this._jsPlumb.types.splice(idx, 1);
          } else this._jsPlumb.types.push(t[i]);
        }
        _applyTypes(this, params, doNotRepaint);
      }
    },
    applyType: function (t, doNotRepaint) {
      this.setPaintStyle(t.paintStyle, doNotRepaint);
      this.setHoverPaintStyle(t.hoverPaintStyle, doNotRepaint);
      if (t.parameters)
        for (var i in t.parameters) this.setParameter(i, t.parameters[i]);
      this._jsPlumb.paintStyleInUse = this.getPaintStyle();
    },
    setPaintStyle: function (style, doNotRepaint) {
      this._jsPlumb.paintStyle = style;
      this._jsPlumb.paintStyleInUse = this._jsPlumb.paintStyle;
      _updateHoverStyle(this);
      if (!doNotRepaint) this.repaint();
    },
    getPaintStyle: function () {
      return this._jsPlumb.paintStyle;
    },
    setHoverPaintStyle: function (style, doNotRepaint) {
      this._jsPlumb.hoverPaintStyle = style;
      _updateHoverStyle(this);
      if (!doNotRepaint) this.repaint();
    },
    getHoverPaintStyle: function () {
      return this._jsPlumb.hoverPaintStyle;
    },
    destroy: function (force) {
      if (force || this.typeId == null) {
        this.cleanupListeners();
        this.clone = null;
        this._jsPlumb = null;
      }
    },
    isHover: function () {
      return this._jsPlumb.hover;
    },
    setHover: function (hover, ignoreAttachedElements, timestamp) {
      if (
        this._jsPlumb &&
        !this._jsPlumb.instance.currentlyDragging &&
        !this._jsPlumb.instance.isHoverSuspended()
      ) {
        this._jsPlumb.hover = hover;
        var method = hover ? "addClass" : "removeClass";
        if (this.canvas != null) {
          if (this._jsPlumb.instance.hoverClass != null)
            this._jsPlumb.instance[method](
              this.canvas,
              this._jsPlumb.instance.hoverClass
            );
          if (this._jsPlumb.hoverClass != null)
            this._jsPlumb.instance[method](
              this.canvas,
              this._jsPlumb.hoverClass
            );
        }
        if (this._jsPlumb.hoverPaintStyle != null) {
          this._jsPlumb.paintStyleInUse = hover
            ? this._jsPlumb.hoverPaintStyle
            : this._jsPlumb.paintStyle;
          if (!this._jsPlumb.instance.isSuspendDrawing()) {
            timestamp = timestamp || jsPlumbUtil.uuid();
            this.repaint({ timestamp, recalc: false });
          }
        }
        if (this.getAttachedElements && !ignoreAttachedElements)
          _updateAttachedElements(this, hover, jsPlumbUtil.uuid(), this);
      }
    },
  });
  var _jsPlumbInstanceIndex = 0;
  var getInstanceIndex = function () {
    var i = _jsPlumbInstanceIndex + 1;
    _jsPlumbInstanceIndex++;
    return i;
  };
  var jsPlumbInstance$jscomp$0 = (root.jsPlumbInstance = function (_defaults) {
    this.version = "2.15.4";
    this.Defaults = {
      Anchor: "Bottom",
      Anchors: [null, null],
      ConnectionsDetachable: true,
      ConnectionOverlays: [],
      Connector: "Bezier",
      Container: null,
      DoNotThrowErrors: false,
      DragOptions: {},
      DropOptions: {},
      Endpoint: "Dot",
      EndpointOverlays: [],
      Endpoints: [null, null],
      EndpointStyle: { fill: "#456" },
      EndpointStyles: [null, null],
      EndpointHoverStyle: null,
      EndpointHoverStyles: [null, null],
      HoverPaintStyle: null,
      LabelStyle: { color: "black" },
      ListStyle: {},
      LogEnabled: false,
      Overlays: [],
      MaxConnections: 1,
      PaintStyle: { "stroke-width": 4, stroke: "#456" },
      ReattachConnections: false,
      RenderMode: "svg",
      Scope: "jsPlumb_DefaultScope",
    };
    if (_defaults) jsPlumb.extend(this.Defaults, _defaults);
    this.logEnabled = this.Defaults.LogEnabled;
    this._connectionTypes = {};
    this._endpointTypes = {};
    _ju.EventGenerator.apply(this);
    var _currentInstance = this;
    var _instanceIndex = getInstanceIndex();
    var _bb = _currentInstance.bind;
    var _initialDefaults = {};
    var _zoom = 1;
    var _info = function (el) {
      if (el == null) return null;
      else if (el.nodeType === 3 || el.nodeType === 8)
        return { el, text: true };
      else {
        var _el = _currentInstance.getElement(el);
        return {
          el: _el,
          id: _ju.isString(el) && _el == null ? el : _getId(_el),
        };
      }
    };
    this.getInstanceIndex = function () {
      return _instanceIndex;
    };
    this.setZoom = function (z, repaintEverything) {
      _zoom = z;
      _currentInstance.fire("zoom", _zoom);
      if (repaintEverything) _currentInstance.repaintEverything();
      return true;
    };
    this.getZoom = function () {
      return _zoom;
    };
    for (var i$jscomp$1 in this.Defaults)
      _initialDefaults[i$jscomp$1] = this.Defaults[i$jscomp$1];
    var _container;
    var _containerDelegations = [];
    this.unbindContainer = function () {
      if (_container != null && _containerDelegations.length > 0)
        for (var i = 0; i < _containerDelegations.length; i++)
          _currentInstance.off(
            _container,
            _containerDelegations[i][0],
            _containerDelegations[i][1]
          );
    };
    this.setContainer = function (c) {
      this.unbindContainer();
      c = this.getElement(c);
      this.select().each(function (conn) {
        conn.moveParent(c);
      });
      this.selectEndpoints().each(function (ep) {
        ep.moveParent(c);
      });
      var previousContainer = _container;
      _container = c;
      _containerDelegations.length = 0;
      var eventAliases = {
        endpointclick: "endpointClick",
        endpointdblclick: "endpointDblClick",
      };
      var _oneDelegateHandler = function (id, e, componentType) {
        var t = e.srcElement || e.target;
        var jp =
          (t && t.parentNode ? t.parentNode._jsPlumb : null) ||
          (t ? t._jsPlumb : null) ||
          (t && t.parentNode && t.parentNode.parentNode
            ? t.parentNode.parentNode._jsPlumb
            : null);
        if (jp) {
          jp.fire(id, jp, e);
          var alias = componentType
            ? eventAliases[componentType + id] || id
            : id;
          _currentInstance.fire(alias, jp.component || jp, e);
        }
      };
      var _addOneDelegate = function (eventId, selector, fn) {
        _containerDelegations.push([eventId, fn]);
        _currentInstance.on(_container, eventId, selector, fn);
      };
      var _oneDelegate = function (id) {
        _addOneDelegate(id, ".jtk-connector", function (e) {
          _oneDelegateHandler(id, e);
        });
        _addOneDelegate(id, ".jtk-endpoint", function (e) {
          _oneDelegateHandler(id, e, "endpoint");
        });
        _addOneDelegate(id, ".jtk-overlay", function (e) {
          _oneDelegateHandler(id, e);
        });
      };
      for (var i = 0; i < events.length; i++) _oneDelegate(events[i]);
      for (var elId in managedElements) {
        var el = managedElements[elId].el;
        if (el.parentNode === previousContainer) {
          previousContainer.removeChild(el);
          _container.appendChild(el);
        }
      }
    };
    this.getContainer = function () {
      return _container;
    };
    this.bind = function (event, fn) {
      if ("ready" === event && initialized) fn();
      else _bb.apply(_currentInstance, [event, fn]);
    };
    _currentInstance.importDefaults = function (d) {
      for (var i in d) _currentInstance.Defaults[i] = d[i];
      if (d.Container) _currentInstance.setContainer(d.Container);
      return _currentInstance;
    };
    _currentInstance.restoreDefaults = function () {
      _currentInstance.Defaults = jsPlumb.extend({}, _initialDefaults);
      return _currentInstance;
    };
    var log = null;
    var initialized = false;
    var connections = [];
    var endpointsByElement = {};
    var endpointsByUUID = {};
    var managedElements = {};
    var offsets = {};
    var offsetTimestamps = {};
    var draggableStates = {};
    var connectionBeingDragged = false;
    var sizes = [];
    var _suspendDrawing = false;
    var _suspendedAt = null;
    var DEFAULT_SCOPE = this.Defaults.Scope;
    var _curIdStamp = 1;
    var _idstamp = function () {
      return "" + _curIdStamp++;
    };
    var _appendElement = function (el, parent) {
      if (_container) _container.appendChild(el);
      else if (!parent) this.appendToRoot(el);
      else this.getElement(parent).appendChild(el);
    }.bind(this);
    var _draw = function (element, ui, timestamp, clearEdits) {
      var drawResult = { c: [], e: [] };
      if (!_suspendDrawing) {
        element = _currentInstance.getElement(element);
        if (element != null) {
          var id = _getId(element);
          var repaintEls = element.querySelectorAll(".jtk-managed");
          if (timestamp == null) timestamp = jsPlumbUtil.uuid();
          var o = _updateOffset({
            elId: id,
            offset: ui,
            recalc: false,
            timestamp,
          });
          for (var i = 0; i < repaintEls.length; i++)
            _updateOffset({
              elId: repaintEls[i].getAttribute("id"),
              recalc: true,
              timestamp,
            });
          var d2 = _currentInstance.router.redraw(
            id,
            ui,
            timestamp,
            null,
            clearEdits
          );
          Array.prototype.push.apply(drawResult.c, d2.c);
          Array.prototype.push.apply(drawResult.e, d2.e);
          if (repaintEls)
            for (var j = 0; j < repaintEls.length; j++) {
              d2 = _currentInstance.router.redraw(
                repaintEls[j].getAttribute("id"),
                null,
                timestamp,
                null,
                clearEdits,
                true
              );
              Array.prototype.push.apply(drawResult.c, d2.c);
              Array.prototype.push.apply(drawResult.e, d2.e);
            }
        }
      }
      return drawResult;
    };
    var _getEndpoint = function (uuid) {
      return endpointsByUUID[uuid];
    };
    var _scopeMatch = function (e1, e2) {
      var s1 = e1.scope.split(/\s/);
      var s2 = e2.scope.split(/\s/);
      for (var i = 0; i < s1.length; i++)
        for (var j = 0; j < s2.length; j++) if (s2[j] === s1[i]) return true;
      return false;
    };
    var _mergeOverrides = function (def, values) {
      var m = jsPlumb.extend({}, def);
      for (var i in values) if (values[i]) m[i] = values[i];
      return m;
    };
    var _prepareConnectionParams = function (params$jscomp$0, referenceParams) {
      var _p = jsPlumb.extend({}, params$jscomp$0);
      if (referenceParams) jsPlumb.extend(_p, referenceParams);
      if (_p.source)
        if (_p.source.endpoint) _p.sourceEndpoint = _p.source;
        else _p.source = _currentInstance.getElement(_p.source);
      if (_p.target)
        if (_p.target.endpoint) _p.targetEndpoint = _p.target;
        else _p.target = _currentInstance.getElement(_p.target);
      if (params$jscomp$0.uuids) {
        _p.sourceEndpoint = _getEndpoint(params$jscomp$0.uuids[0]);
        _p.targetEndpoint = _getEndpoint(params$jscomp$0.uuids[1]);
      }
      if (_p.sourceEndpoint && _p.sourceEndpoint.isFull()) {
        _ju.log(
          _currentInstance,
          "could not add connection; source endpoint is full"
        );
        return;
      }
      if (_p.targetEndpoint && _p.targetEndpoint.isFull()) {
        _ju.log(
          _currentInstance,
          "could not add connection; target endpoint is full"
        );
        return;
      }
      if (!_p.type && _p.sourceEndpoint)
        _p.type = _p.sourceEndpoint.connectionType;
      if (_p.sourceEndpoint && _p.sourceEndpoint.connectorOverlays) {
        _p.overlays = _p.overlays || [];
        var i = 0;
        for (var j = _p.sourceEndpoint.connectorOverlays.length; i < j; i++)
          _p.overlays.push(_p.sourceEndpoint.connectorOverlays[i]);
      }
      if (_p.sourceEndpoint && _p.sourceEndpoint.scope)
        _p.scope = _p.sourceEndpoint.scope;
      if (
        !_p["pointer-events"] &&
        _p.sourceEndpoint &&
        _p.sourceEndpoint.connectorPointerEvents
      )
        _p["pointer-events"] = _p.sourceEndpoint.connectorPointerEvents;
      var _addEndpoint = function (el, def, idx) {
        var params = _mergeOverrides(def, {
          anchor: _p.anchors ? _p.anchors[idx] : _p.anchor,
          endpoint: _p.endpoints ? _p.endpoints[idx] : _p.endpoint,
          paintStyle: _p.endpointStyles
            ? _p.endpointStyles[idx]
            : _p.endpointStyle,
          hoverPaintStyle: _p.endpointHoverStyles
            ? _p.endpointHoverStyles[idx]
            : _p.endpointHoverStyle,
        });
        return _currentInstance.addEndpoint(el, params);
      };
      var _oneElementDef = function (type, idx, defs, matchType) {
        if (
          _p[type] &&
          !_p[type].endpoint &&
          !_p[type + "Endpoint"] &&
          !_p.newConnection
        ) {
          var tid = _getId(_p[type]);
          var tep = defs[tid];
          tep = tep ? tep[matchType] : null;
          if (tep) {
            if (!tep.enabled) return false;
            var epDef = jsPlumb.extend({}, tep.def);
            delete epDef.label;
            var newEndpoint =
              tep.endpoint != null && tep.endpoint._jsPlumb
                ? tep.endpoint
                : _addEndpoint(_p[type], epDef, idx);
            if (newEndpoint.isFull()) return false;
            _p[type + "Endpoint"] = newEndpoint;
            if (!_p.scope && epDef.scope) _p.scope = epDef.scope;
            if (tep.uniqueEndpoint)
              if (!tep.endpoint) {
                tep.endpoint = newEndpoint;
                newEndpoint.setDeleteOnEmpty(false);
              } else newEndpoint.finalEndpoint = tep.endpoint;
            else newEndpoint.setDeleteOnEmpty(true);
            if (idx === 0 && tep.def.connectorOverlays) {
              _p.overlays = _p.overlays || [];
              Array.prototype.push.apply(
                _p.overlays,
                tep.def.connectorOverlays
              );
            }
          }
        }
      };
      if (
        _oneElementDef(
          "source",
          0,
          this.sourceEndpointDefinitions,
          _p.type || "default"
        ) === false
      )
        return;
      if (
        _oneElementDef(
          "target",
          1,
          this.targetEndpointDefinitions,
          _p.type || "default"
        ) === false
      )
        return;
      if (_p.sourceEndpoint && _p.targetEndpoint)
        if (!_scopeMatch(_p.sourceEndpoint, _p.targetEndpoint)) _p = null;
      return _p;
    }.bind(_currentInstance);
    var _newConnection = function (params) {
      var connectionFunc =
        _currentInstance.Defaults.ConnectionType ||
        _currentInstance.getDefaultConnectionType();
      params._jsPlumb = _currentInstance;
      params.newConnection = _newConnection;
      params.newEndpoint = _newEndpoint;
      params.endpointsByUUID = endpointsByUUID;
      params.endpointsByElement = endpointsByElement;
      params.finaliseConnection = _finaliseConnection;
      params.id = "con_" + _idstamp();
      var con = new connectionFunc(params);
      if (con.isDetachable()) {
        con.endpoints[0].initDraggable("_jsPlumbSource");
        con.endpoints[1].initDraggable("_jsPlumbTarget");
      }
      return con;
    };
    var _finaliseConnection = (_currentInstance.finaliseConnection = function (
      jpc,
      params,
      originalEvent,
      doInformAnchorManager
    ) {
      params = params || {};
      if (!jpc.suspendedEndpoint) connections.push(jpc);
      jpc.pending = null;
      jpc.endpoints[0].isTemporarySource = false;
      if (doInformAnchorManager !== false)
        _currentInstance.router.newConnection(jpc);
      _draw(jpc.source);
      if (!params.doNotFireConnectionEvent && params.fireEvent !== false) {
        var eventArgs = {
          connection: jpc,
          source: jpc.source,
          target: jpc.target,
          sourceId: jpc.sourceId,
          targetId: jpc.targetId,
          sourceEndpoint: jpc.endpoints[0],
          targetEndpoint: jpc.endpoints[1],
        };
        _currentInstance.fire("connection", eventArgs, originalEvent);
      }
    });
    var _newEndpoint = function (params, id) {
      var endpointFunc =
        _currentInstance.Defaults.EndpointType || jsPlumb.Endpoint;
      var _p = jsPlumb.extend({}, params);
      _p._jsPlumb = _currentInstance;
      _p.newConnection = _newConnection;
      _p.newEndpoint = _newEndpoint;
      _p.endpointsByUUID = endpointsByUUID;
      _p.endpointsByElement = endpointsByElement;
      _p.fireDetachEvent = fireDetachEvent;
      _p.elementId = id || _getId(_p.source);
      var ep = new endpointFunc(_p);
      ep.id = "ep_" + _idstamp();
      _manage(_p.elementId, _p.source);
      if (!jsPlumb.headless)
        _currentInstance.getDragManager().endpointAdded(_p.source, id);
      return ep;
    };
    var _operation = function (elId, func, endpointFunc) {
      var endpoints = endpointsByElement[elId];
      if (endpoints && endpoints.length) {
        var i = 0;
        for (var ii = endpoints.length; i < ii; i++) {
          var j = 0;
          for (var jj = endpoints[i].connections.length; j < jj; j++) {
            var retVal = func(endpoints[i].connections[j]);
            if (retVal) return;
          }
          if (endpointFunc) endpointFunc(endpoints[i]);
        }
      }
    };
    var _setDraggable = function (element, draggable) {
      return jsPlumb.each(element, function (el) {
        if (_currentInstance.isDragSupported(el)) {
          draggableStates[_currentInstance.getAttribute(el, "id")] = draggable;
          _currentInstance.setElementDraggable(el, draggable);
        }
      });
    };
    var _setVisible = function (el, state, alsoChangeEndpoints) {
      state = state === "block";
      var endpointFunc = null;
      if (alsoChangeEndpoints)
        endpointFunc = function (ep) {
          ep.setVisible(state, true, true);
        };
      var info = _info(el);
      _operation(
        info.id,
        function (jpc) {
          if (state && alsoChangeEndpoints) {
            var oidx = jpc.sourceId === info.id ? 1 : 0;
            if (jpc.endpoints[oidx].isVisible()) jpc.setVisible(true);
          } else jpc.setVisible(state);
        },
        endpointFunc
      );
    };
    var _toggleVisible = function (elId, changeEndpoints) {
      var endpointFunc = null;
      if (changeEndpoints)
        endpointFunc = function (ep) {
          var state = ep.isVisible();
          ep.setVisible(!state);
        };
      _operation(
        elId,
        function (jpc) {
          var state = jpc.isVisible();
          jpc.setVisible(!state);
        },
        endpointFunc
      );
    };
    var _getCachedData = function (elId) {
      var o = offsets[elId];
      if (!o) return _updateOffset({ elId });
      else return { o, s: sizes[elId] };
    };
    var _getId = function (element, uuid, doNotCreateIfNotFound) {
      if (_ju.isString(element)) return element;
      if (element == null) return null;
      var id = _currentInstance.getAttribute(element, "id");
      if (!id || id === "undefined") {
        if (arguments.length === 2 && arguments[1] !== undefined) id = uuid;
        else if (
          arguments.length === 1 ||
          (arguments.length === 3 && !arguments[2])
        )
          id = "jsPlumb_" + _instanceIndex + "_" + _idstamp();
        if (!doNotCreateIfNotFound)
          _currentInstance.setAttribute(element, "id", id);
      }
      return id;
    };
    this.setConnectionBeingDragged = function (v) {
      connectionBeingDragged = v;
    };
    this.isConnectionBeingDragged = function () {
      return connectionBeingDragged;
    };
    this.getManagedElements = function () {
      return managedElements;
    };
    this.connectorClass = "jtk-connector";
    this.connectorOutlineClass = "jtk-connector-outline";
    this.connectedClass = "jtk-connected";
    this.hoverClass = "jtk-hover";
    this.endpointClass = "jtk-endpoint";
    this.endpointConnectedClass = "jtk-endpoint-connected";
    this.endpointFullClass = "jtk-endpoint-full";
    this.endpointDropAllowedClass = "jtk-endpoint-drop-allowed";
    this.endpointDropForbiddenClass = "jtk-endpoint-drop-forbidden";
    this.overlayClass = "jtk-overlay";
    this.draggingClass = "jtk-dragging";
    this.elementDraggingClass = "jtk-element-dragging";
    this.sourceElementDraggingClass = "jtk-source-element-dragging";
    this.targetElementDraggingClass = "jtk-target-element-dragging";
    this.endpointAnchorClassPrefix = "jtk-endpoint-anchor";
    this.hoverSourceClass = "jtk-source-hover";
    this.hoverTargetClass = "jtk-target-hover";
    this.dragSelectClass = "jtk-drag-select";
    this.Anchors = {};
    this.Connectors = { svg: {} };
    this.Endpoints = { svg: {} };
    this.Overlays = { svg: {} };
    this.ConnectorRenderers = {};
    this.SVG = "svg";
    this.addEndpoint = function (el, params, referenceParams) {
      referenceParams = referenceParams || {};
      var p = jsPlumb.extend({}, referenceParams);
      jsPlumb.extend(p, params);
      p.endpoint = p.endpoint || _currentInstance.Defaults.Endpoint;
      p.paintStyle = p.paintStyle || _currentInstance.Defaults.EndpointStyle;
      var results = [];
      var inputs =
        _ju.isArray(el) || (el.length != null && !_ju.isString(el)) ? el : [el];
      var i = 0;
      for (var j = inputs.length; i < j; i++) {
        p.source = _currentInstance.getElement(inputs[i]);
        _ensureContainer(p.source);
        var id = _getId(p.source);
        var e = _newEndpoint(p, id);
        var myOffset = _manage(id, p.source, null, !_suspendDrawing).info.o;
        _ju.addToList(endpointsByElement, id, e);
        if (!_suspendDrawing)
          e.paint({
            anchorLoc: e.anchor.compute({
              xy: [myOffset.left, myOffset.top],
              wh: sizes[id],
              element: e,
              timestamp: _suspendedAt,
              rotation: this.getRotation(id),
            }),
            timestamp: _suspendedAt,
          });
        results.push(e);
      }
      return results.length === 1 ? results[0] : results;
    };
    this.addEndpoints = function (el, endpoints, referenceParams) {
      var results = [];
      var i = 0;
      for (var j = endpoints.length; i < j; i++) {
        var e = _currentInstance.addEndpoint(el, endpoints[i], referenceParams);
        if (_ju.isArray(e)) Array.prototype.push.apply(results, e);
        else results.push(e);
      }
      return results;
    };
    this.animate = function (el, properties, options) {
      if (!this.animationSupported) return false;
      options = options || {};
      var del = _currentInstance.getElement(el);
      var id = _getId(del);
      var stepFunction = jsPlumb.animEvents.step;
      var completeFunction = jsPlumb.animEvents.complete;
      options[stepFunction] = _ju.wrap(options[stepFunction], function () {
        _currentInstance.revalidate(id);
      });
      options[completeFunction] = _ju.wrap(
        options[completeFunction],
        function () {
          _currentInstance.revalidate(id);
        }
      );
      _currentInstance.doAnimate(del, properties, options);
    };
    this.checkCondition = function (conditionName, args) {
      var l = _currentInstance.getListener(conditionName);
      var r = true;
      if (l && l.length > 0) {
        var values = Array.prototype.slice.call(arguments, 1);
        try {
          var i = 0;
          for (var j = l.length; i < j; i++) r = r && l[i].apply(l[i], values);
        } catch (e) {
          _ju.log(
            _currentInstance,
            "cannot check condition [" + conditionName + "]" + e
          );
        }
      }
      return r;
    };
    this.connect = function (params, referenceParams) {
      var _p = _prepareConnectionParams(params, referenceParams);
      if (_p) {
        if (_p.source == null && _p.sourceEndpoint == null) {
          _ju.log("Cannot establish connection - source does not exist");
          return;
        }
        if (_p.target == null && _p.targetEndpoint == null) {
          _ju.log("Cannot establish connection - target does not exist");
          return;
        }
        _ensureContainer(_p.source);
        var jpc = _newConnection(_p);
        _finaliseConnection(jpc, _p);
      }
      return jpc;
    };
    var stTypes = [
      { el: "source", elId: "sourceId", epDefs: "sourceEndpointDefinitions" },
      { el: "target", elId: "targetId", epDefs: "targetEndpointDefinitions" },
    ];
    var _set = function (c, el, idx, doNotRepaint) {
      var _st = stTypes[idx];
      var cId = c[_st.elId];
      var cEl = c[_st.el];
      var oldEndpoint = c.endpoints[idx];
      var evtParams = {
        index: idx,
        originalSourceId: idx === 0 ? cId : c.sourceId,
        newSourceId: c.sourceId,
        originalTargetId: idx === 1 ? cId : c.targetId,
        newTargetId: c.targetId,
        connection: c,
      };
      if (el.constructor === jsPlumb.Endpoint) {
        var ep = el;
        ep.addConnection(c);
        el = ep.element;
      } else {
        var sid = _getId(el);
        var sep = this[_st.epDefs][sid];
        if (sid === c[_st.elId]) ep = null;
        else if (sep)
          for (var t in sep) {
            if (!sep[t].enabled) return;
            ep =
              sep[t].endpoint != null && sep[t].endpoint._jsPlumb
                ? sep[t].endpoint
                : this.addEndpoint(el, sep[t].def);
            if (sep[t].uniqueEndpoint) sep[t].endpoint = ep;
            ep.addConnection(c);
          }
        else ep = c.makeEndpoint(idx === 0, el, sid);
      }
      if (ep != null) {
        oldEndpoint.detachFromConnection(c);
        c.endpoints[idx] = ep;
        c[_st.el] = ep.element;
        c[_st.elId] = ep.elementId;
        evtParams[idx === 0 ? "newSourceId" : "newTargetId"] = ep.elementId;
        fireMoveEvent(evtParams);
        if (!doNotRepaint) c.repaint();
      }
      evtParams.element = el;
      return evtParams;
    }.bind(this);
    this.setSource = function (connection, el, doNotRepaint) {
      var p = _set(connection, el, 0, doNotRepaint);
      this.router.sourceOrTargetChanged(
        p.originalSourceId,
        p.newSourceId,
        connection,
        p.el,
        0
      );
    };
    this.setTarget = function (connection, el, doNotRepaint) {
      var p = _set(connection, el, 1, doNotRepaint);
      this.router.sourceOrTargetChanged(
        p.originalTargetId,
        p.newTargetId,
        connection,
        p.el,
        1
      );
    };
    this.deleteEndpoint = function (
      object,
      dontUpdateHover,
      deleteAttachedObjects
    ) {
      var endpoint =
        typeof object === "string" ? endpointsByUUID[object] : object;
      if (endpoint)
        _currentInstance.deleteObject({
          endpoint,
          dontUpdateHover,
          deleteAttachedObjects,
        });
      return _currentInstance;
    };
    this.deleteEveryEndpoint = function () {
      var _is = _currentInstance.setSuspendDrawing(true);
      for (var id in endpointsByElement) {
        var endpoints = endpointsByElement[id];
        if (endpoints && endpoints.length) {
          var i = 0;
          for (var j = endpoints.length; i < j; i++)
            _currentInstance.deleteEndpoint(endpoints[i], true);
        }
      }
      endpointsByElement = {};
      managedElements = {};
      endpointsByUUID = {};
      offsets = {};
      offsetTimestamps = {};
      _currentInstance.router.reset();
      var dm = _currentInstance.getDragManager();
      if (dm) dm.reset();
      if (!_is) _currentInstance.setSuspendDrawing(false);
      return _currentInstance;
    };
    var fireDetachEvent = function (jpc, doFireEvent, originalEvent) {
      var connType =
        _currentInstance.Defaults.ConnectionType ||
        _currentInstance.getDefaultConnectionType();
      var argIsConnection = jpc.constructor === connType;
      var params = argIsConnection
        ? {
            connection: jpc,
            source: jpc.source,
            target: jpc.target,
            sourceId: jpc.sourceId,
            targetId: jpc.targetId,
            sourceEndpoint: jpc.endpoints[0],
            targetEndpoint: jpc.endpoints[1],
          }
        : jpc;
      if (doFireEvent)
        _currentInstance.fire("connectionDetached", params, originalEvent);
      _currentInstance.fire(
        "internal.connectionDetached",
        params,
        originalEvent
      );
      _currentInstance.router.connectionDetached(params);
    };
    var fireMoveEvent = (_currentInstance.fireMoveEvent = function (
      params,
      evt
    ) {
      _currentInstance.fire("connectionMoved", params, evt);
    });
    this.unregisterEndpoint = function (endpoint) {
      if (endpoint._jsPlumb.uuid)
        endpointsByUUID[endpoint._jsPlumb.uuid] = null;
      _currentInstance.router.deleteEndpoint(endpoint);
      for (var e in endpointsByElement) {
        var endpoints = endpointsByElement[e];
        if (endpoints) {
          var newEndpoints = [];
          var i = 0;
          for (var j = endpoints.length; i < j; i++)
            if (endpoints[i] !== endpoint) newEndpoints.push(endpoints[i]);
          endpointsByElement[e] = newEndpoints;
        }
        if (endpointsByElement[e].length < 1) delete endpointsByElement[e];
      }
    };
    var IS_DETACH_ALLOWED = "isDetachAllowed";
    var BEFORE_DETACH = "beforeDetach";
    var CHECK_CONDITION = "checkCondition";
    this.deleteConnection = function (connection, params) {
      if (connection != null) {
        params = params || {};
        if (
          params.force ||
          _ju.functionChain(true, false, [
            [connection.endpoints[0], IS_DETACH_ALLOWED, [connection]],
            [connection.endpoints[1], IS_DETACH_ALLOWED, [connection]],
            [connection, IS_DETACH_ALLOWED, [connection]],
            [_currentInstance, CHECK_CONDITION, [BEFORE_DETACH, connection]],
          ])
        ) {
          connection.setHover(false);
          fireDetachEvent(
            connection,
            !connection.pending && params.fireEvent !== false,
            params.originalEvent
          );
          connection.endpoints[0].detachFromConnection(connection);
          connection.endpoints[1].detachFromConnection(connection);
          _ju.removeWithFunction(connections, function (_c) {
            return connection.id === _c.id;
          });
          connection.cleanup();
          connection.destroy();
          return true;
        }
      }
      return false;
    };
    this.deleteEveryConnection = function (params) {
      params = params || {};
      var count = connections.length;
      var deletedCount = 0;
      _currentInstance.batch(function () {
        for (var i = 0; i < count; i++)
          deletedCount += _currentInstance.deleteConnection(
            connections[0],
            params
          )
            ? 1
            : 0;
      });
      return deletedCount;
    };
    this.deleteConnectionsForElement = function (el, params) {
      params = params || {};
      el = _currentInstance.getElement(el);
      var id = _getId(el);
      var endpoints = endpointsByElement[id];
      if (endpoints && endpoints.length) {
        var i = 0;
        for (var j = endpoints.length; i < j; i++)
          endpoints[i].deleteEveryConnection(params);
      }
      return _currentInstance;
    };
    this.deleteObject = function (params) {
      var result = {
        endpoints: {},
        connections: {},
        endpointCount: 0,
        connectionCount: 0,
      };
      var deleteAttachedObjects = params.deleteAttachedObjects !== false;
      var unravelConnection = function (connection) {
        if (connection != null && result.connections[connection.id] == null) {
          if (!params.dontUpdateHover && connection._jsPlumb != null)
            connection.setHover(false);
          result.connections[connection.id] = connection;
          result.connectionCount++;
        }
      };
      var unravelEndpoint = function (endpoint) {
        if (endpoint != null && result.endpoints[endpoint.id] == null) {
          if (!params.dontUpdateHover && endpoint._jsPlumb != null)
            endpoint.setHover(false);
          result.endpoints[endpoint.id] = endpoint;
          result.endpointCount++;
          if (deleteAttachedObjects)
            for (var i = 0; i < endpoint.connections.length; i++) {
              var c = endpoint.connections[i];
              unravelConnection(c);
            }
        }
      };
      if (params.connection) unravelConnection(params.connection);
      else unravelEndpoint(params.endpoint);
      for (var i$jscomp$0 in result.connections) {
        var c$jscomp$0 = result.connections[i$jscomp$0];
        if (c$jscomp$0._jsPlumb) {
          _ju.removeWithFunction(connections, function (_c) {
            return c$jscomp$0.id === _c.id;
          });
          fireDetachEvent(
            c$jscomp$0,
            params.fireEvent === false ? false : !c$jscomp$0.pending,
            params.originalEvent
          );
          var doNotCleanup =
            params.deleteAttachedObjects == null
              ? null
              : !params.deleteAttachedObjects;
          c$jscomp$0.endpoints[0].detachFromConnection(
            c$jscomp$0,
            null,
            doNotCleanup
          );
          c$jscomp$0.endpoints[1].detachFromConnection(
            c$jscomp$0,
            null,
            doNotCleanup
          );
          c$jscomp$0.cleanup(true);
          c$jscomp$0.destroy(true);
        }
      }
      for (var j in result.endpoints) {
        var e = result.endpoints[j];
        if (e._jsPlumb) {
          _currentInstance.unregisterEndpoint(e);
          e.cleanup(true);
          e.destroy(true);
        }
      }
      return result;
    };
    var _setOperation = function (list, func, args, selector) {
      var i = 0;
      for (var j = list.length; i < j; i++) list[i][func].apply(list[i], args);
      return selector(list);
    };
    var _getOperation = function (list, func, args) {
      var out = [];
      var i = 0;
      for (var j = list.length; i < j; i++)
        out.push([list[i][func].apply(list[i], args), list[i]]);
      return out;
    };
    var setter = function (list, func, selector) {
      return function () {
        return _setOperation(list, func, arguments, selector);
      };
    };
    var getter = function (list, func) {
      return function () {
        return _getOperation(list, func, arguments);
      };
    };
    var prepareList = function (input, doNotGetIds) {
      var r = [];
      if (input)
        if (typeof input === "string") {
          if (input === "*") return input;
          r.push(input);
        } else if (doNotGetIds) r = input;
        else if (input.length) {
          var i = 0;
          for (var j = input.length; i < j; i++) r.push(_info(input[i]).id);
        } else r.push(_info(input).id);
      return r;
    };
    var filterList = function (list, value, missingIsFalse) {
      if (list === "*") return true;
      return list.length > 0 ? list.indexOf(value) !== -1 : !missingIsFalse;
    };
    this.getConnections = function (options, flat) {
      if (!options) options = {};
      else if (options.constructor === String) options = { scope: options };
      var scope$jscomp$0 = options.scope || _currentInstance.getDefaultScope();
      var scopes = prepareList(scope$jscomp$0, true);
      var sources = prepareList(options.source);
      var targets = prepareList(options.target);
      var results = !flat && scopes.length > 1 ? {} : [];
      var _addOne = function (scope, obj) {
        if (!flat && scopes.length > 1) {
          var ss = results[scope];
          if (ss == null) ss = results[scope] = [];
          ss.push(obj);
        } else results.push(obj);
      };
      var j = 0;
      for (var jj = connections.length; j < jj; j++) {
        var c = connections[j];
        var sourceId =
          c.proxies && c.proxies[0]
            ? c.proxies[0].originalEp.elementId
            : c.sourceId;
        var targetId =
          c.proxies && c.proxies[1]
            ? c.proxies[1].originalEp.elementId
            : c.targetId;
        if (
          filterList(scopes, c.scope) &&
          filterList(sources, sourceId) &&
          filterList(targets, targetId)
        )
          _addOne(c.scope, c);
      }
      return results;
    };
    var _curryEach = function (list, executor) {
      return function (f) {
        var i = 0;
        for (var ii = list.length; i < ii; i++) f(list[i]);
        return executor(list);
      };
    };
    var _curryGet = function (list) {
      return function (idx) {
        return list[idx];
      };
    };
    var _makeCommonSelectHandler = function (list, executor) {
      var out = {
        length: list.length,
        each: _curryEach(list, executor),
        get: _curryGet(list),
      };
      var setters = [
        "setHover",
        "removeAllOverlays",
        "setLabel",
        "addClass",
        "addOverlay",
        "removeOverlay",
        "removeOverlays",
        "showOverlay",
        "hideOverlay",
        "showOverlays",
        "hideOverlays",
        "setPaintStyle",
        "setHoverPaintStyle",
        "setSuspendEvents",
        "setParameter",
        "setParameters",
        "setVisible",
        "repaint",
        "addType",
        "toggleType",
        "removeType",
        "removeClass",
        "setType",
        "bind",
        "unbind",
      ];
      var getters = [
        "getLabel",
        "getOverlay",
        "isHover",
        "getParameter",
        "getParameters",
        "getPaintStyle",
        "getHoverPaintStyle",
        "isVisible",
        "hasType",
        "getType",
        "isSuspendEvents",
      ];
      var i;
      var ii;
      for (i = 0, ii = setters.length; i < ii; i++)
        out[setters[i]] = setter(list, setters[i], executor);
      for (i = 0, ii = getters.length; i < ii; i++)
        out[getters[i]] = getter(list, getters[i]);
      return out;
    };
    var _makeConnectionSelectHandler = function (list) {
      var common = _makeCommonSelectHandler(list, _makeConnectionSelectHandler);
      return jsPlumb.extend(common, {
        setDetachable: setter(
          list,
          "setDetachable",
          _makeConnectionSelectHandler
        ),
        setReattach: setter(list, "setReattach", _makeConnectionSelectHandler),
        setConnector: setter(
          list,
          "setConnector",
          _makeConnectionSelectHandler
        ),
        delete: function () {
          var i = 0;
          for (var ii = list.length; i < ii; i++)
            _currentInstance.deleteConnection(list[i]);
        },
        isDetachable: getter(list, "isDetachable"),
        isReattach: getter(list, "isReattach"),
      });
    };
    var _makeEndpointSelectHandler = function (list) {
      var common = _makeCommonSelectHandler(list, _makeEndpointSelectHandler);
      return jsPlumb.extend(common, {
        setEnabled: setter(list, "setEnabled", _makeEndpointSelectHandler),
        setAnchor: setter(list, "setAnchor", _makeEndpointSelectHandler),
        isEnabled: getter(list, "isEnabled"),
        deleteEveryConnection: function () {
          var i = 0;
          for (var ii = list.length; i < ii; i++)
            list[i].deleteEveryConnection();
        },
        delete: function () {
          var i = 0;
          for (var ii = list.length; i < ii; i++)
            _currentInstance.deleteEndpoint(list[i]);
        },
      });
    };
    this.select = function (params) {
      params = params || {};
      params.scope = params.scope || "*";
      return _makeConnectionSelectHandler(
        params.connections || _currentInstance.getConnections(params, true)
      );
    };
    this.selectEndpoints = function (params) {
      params = params || {};
      params.scope = params.scope || "*";
      var noElementFilters =
        !params.element && !params.source && !params.target;
      var elements = noElementFilters ? "*" : prepareList(params.element);
      var sources = noElementFilters ? "*" : prepareList(params.source);
      var targets = noElementFilters ? "*" : prepareList(params.target);
      var scopes = prepareList(params.scope, true);
      var ep = [];
      for (var el in endpointsByElement) {
        var either = filterList(elements, el, true);
        var source = filterList(sources, el, true);
        var sourceMatchExact = sources !== "*";
        var target = filterList(targets, el, true);
        var targetMatchExact = targets !== "*";
        if (either || source || target) {
          var i = 0;
          var ii = endpointsByElement[el].length;
          inner: for (; i < ii; i++) {
            var _ep = endpointsByElement[el][i];
            if (filterList(scopes, _ep.scope, true)) {
              var noMatchSource =
                sourceMatchExact && sources.length > 0 && !_ep.isSource;
              var noMatchTarget =
                targetMatchExact && targets.length > 0 && !_ep.isTarget;
              if (noMatchSource || noMatchTarget) continue inner;
              ep.push(_ep);
            }
          }
        }
      }
      return _makeEndpointSelectHandler(ep);
    };
    this.getAllConnections = function () {
      return connections;
    };
    this.getDefaultScope = function () {
      return DEFAULT_SCOPE;
    };
    this.getEndpoint = _getEndpoint;
    this.getEndpoints = function (el) {
      return endpointsByElement[_info(el).id] || [];
    };
    this.getDefaultEndpointType = function () {
      return jsPlumb.Endpoint;
    };
    this.getDefaultConnectionType = function () {
      return jsPlumb.Connection;
    };
    this.getId = _getId;
    this.draw = _draw;
    this.info = _info;
    this.appendElement = _appendElement;
    var _hoverSuspended = false;
    this.isHoverSuspended = function () {
      return _hoverSuspended;
    };
    this.setHoverSuspended = function (s) {
      _hoverSuspended = s;
    };
    this.hide = function (el, changeEndpoints) {
      _setVisible(el, "none", changeEndpoints);
      return _currentInstance;
    };
    this.idstamp = _idstamp;
    var _ensureContainer = function (candidate) {
      if (!_container && candidate) {
        var can = _currentInstance.getElement(candidate);
        if (can.offsetParent) _currentInstance.setContainer(can.offsetParent);
      }
    };
    var _getContainerFromDefaults = function () {
      if (_currentInstance.Defaults.Container)
        _currentInstance.setContainer(_currentInstance.Defaults.Container);
    };
    var _manage = (_currentInstance.manage = function (
      id,
      element,
      _transient,
      _recalc
    ) {
      if (!managedElements[id]) {
        managedElements[id] = {
          el: element,
          endpoints: [],
          connections: [],
          rotation: 0,
        };
        managedElements[id].info = _updateOffset({
          elId: id,
          timestamp: _suspendedAt,
        });
        _currentInstance.addClass(element, "jtk-managed");
        if (!_transient)
          _currentInstance.fire("manageElement", {
            id,
            info: managedElements[id].info,
            el: element,
          });
      } else if (_recalc)
        managedElements[id].info = _updateOffset({
          elId: id,
          timestamp: _suspendedAt,
          recalc: true,
        });
      return managedElements[id];
    });
    this.unmanage = function (id) {
      if (managedElements[id]) {
        var el = managedElements[id].el;
        _currentInstance.removeClass(el, "jtk-managed");
        delete managedElements[id];
        _currentInstance.fire("unmanageElement", { id, el });
      }
    };
    this.rotate = function (elId, amountInDegrees, doNotRedraw) {
      if (managedElements[elId]) {
        managedElements[elId].rotation = amountInDegrees;
        managedElements[elId].el.style.transform =
          "rotate(" + amountInDegrees + "deg)";
        managedElements[elId].el.style.transformOrigin = "center center";
        if (doNotRedraw !== true) return this.revalidate(elId);
      }
      return { c: [], e: [] };
    };
    this.getRotation = function (elementId) {
      return managedElements[elementId]
        ? managedElements[elementId].rotation || 0
        : 0;
    };
    var _updateOffset = function (params) {
      var timestamp = params.timestamp;
      var recalc = params.recalc;
      var offset = params.offset;
      var elId = params.elId;
      if (_suspendDrawing && !timestamp) timestamp = _suspendedAt;
      if (!recalc)
        if (timestamp && timestamp === offsetTimestamps[elId])
          return { o: params.offset || offsets[elId], s: sizes[elId] };
      if (recalc || (!offset && offsets[elId] == null)) {
        var s = managedElements[elId] ? managedElements[elId].el : null;
        if (s != null) {
          sizes[elId] = _currentInstance.getSize(s);
          offsets[elId] = _currentInstance.getOffset(s);
          offsetTimestamps[elId] = timestamp;
        }
      } else {
        offsets[elId] = offset || offsets[elId];
        if (sizes[elId] == null) {
          s = managedElements[elId].el;
          if (s != null) sizes[elId] = _currentInstance.getSize(s);
        }
        offsetTimestamps[elId] = timestamp;
      }
      if (offsets[elId] && !offsets[elId].right) {
        offsets[elId].right = offsets[elId].left + sizes[elId][0];
        offsets[elId].bottom = offsets[elId].top + sizes[elId][1];
        offsets[elId].width = sizes[elId][0];
        offsets[elId].height = sizes[elId][1];
        offsets[elId].centerx = offsets[elId].left + offsets[elId].width / 2;
        offsets[elId].centery = offsets[elId].top + offsets[elId].height / 2;
      }
      return { o: offsets[elId], s: sizes[elId] };
    };
    this.updateOffset = _updateOffset;
    this.init = function () {
      if (!initialized) {
        _getContainerFromDefaults();
        _currentInstance.router = new root.jsPlumb.DefaultRouter(
          _currentInstance
        );
        _currentInstance.anchorManager = _currentInstance.router.anchorManager;
        initialized = true;
        _currentInstance.fire("ready", _currentInstance);
      }
    }.bind(this);
    this.log = log;
    this.jsPlumbUIComponent = jsPlumbUIComponent;
    this.makeAnchor = function () {
      var _a = function (t, p) {
        if (root.jsPlumb.Anchors[t]) return new root.jsPlumb.Anchors[t](p);
        if (!_currentInstance.Defaults.DoNotThrowErrors)
          throw { msg: "jsPlumb: unknown anchor type '" + t + "'" };
      };
      if (arguments.length === 0) return null;
      var specimen = arguments[0];
      var elementId = arguments[1];
      var jsPlumbInstance = arguments[2];
      var newAnchor = null;
      if (specimen.compute && specimen.getOrientation) return specimen;
      else if (typeof specimen === "string")
        newAnchor = _a(arguments[0], {
          elementId,
          jsPlumbInstance: _currentInstance,
        });
      else if (_ju.isArray(specimen))
        if (_ju.isArray(specimen[0]) || _ju.isString(specimen[0]))
          if (specimen.length === 2 && _ju.isObject(specimen[1]))
            if (_ju.isString(specimen[0])) {
              var pp = root.jsPlumb.extend(
                { elementId, jsPlumbInstance: _currentInstance },
                specimen[1]
              );
              newAnchor = _a(specimen[0], pp);
            } else {
              pp = root.jsPlumb.extend(
                {
                  elementId,
                  jsPlumbInstance: _currentInstance,
                  anchors: specimen[0],
                },
                specimen[1]
              );
              newAnchor = new root.jsPlumb.DynamicAnchor(pp);
            }
          else
            newAnchor = new jsPlumb.DynamicAnchor({
              anchors: specimen,
              selector: null,
              elementId,
              jsPlumbInstance: _currentInstance,
            });
        else {
          var anchorParams = {
            x: specimen[0],
            y: specimen[1],
            orientation:
              specimen.length >= 4 ? [specimen[2], specimen[3]] : [0, 0],
            offsets: specimen.length >= 6 ? [specimen[4], specimen[5]] : [0, 0],
            elementId,
            jsPlumbInstance: _currentInstance,
            cssClass: specimen.length === 7 ? specimen[6] : null,
          };
          newAnchor = new root.jsPlumb.Anchor(anchorParams);
          newAnchor.clone = function () {
            return new root.jsPlumb.Anchor(anchorParams);
          };
        }
      if (!newAnchor.id) newAnchor.id = "anchor_" + _idstamp();
      return newAnchor;
    };
    this.makeAnchors = function (types, elementId, jsPlumbInstance) {
      var r = [];
      var i = 0;
      for (var ii = types.length; i < ii; i++)
        if (typeof types[i] === "string")
          r.push(
            root.jsPlumb.Anchors[types[i]]({ elementId, jsPlumbInstance })
          );
        else if (_ju.isArray(types[i]))
          r.push(
            _currentInstance.makeAnchor(types[i], elementId, jsPlumbInstance)
          );
      return r;
    };
    this.makeDynamicAnchor = function (anchors, anchorSelector) {
      return new root.jsPlumb.DynamicAnchor({
        anchors,
        selector: anchorSelector,
        elementId: null,
        jsPlumbInstance: _currentInstance,
      });
    };
    this.targetEndpointDefinitions = {};
    this.sourceEndpointDefinitions = {};
    var selectorFilter = function (evt, _el, selector, _instance, negate) {
      var t = evt.target || evt.srcElement;
      var ok = false;
      var sel = _instance.getSelector(_el, selector);
      for (var j = 0; j < sel.length; j++)
        if (sel[j] === t) {
          ok = true;
          break;
        }
      return negate ? !ok : ok;
    };
    var _makeElementDropHandler = function (
      elInfo,
      p,
      dropOptions,
      isSource,
      isTarget
    ) {
      var proxyComponent = new jsPlumbUIComponent(p);
      var _drop = p._jsPlumb.EndpointDropHandler({
        jsPlumb: _currentInstance,
        enabled: function () {
          return elInfo.def.enabled;
        },
        isFull: function () {
          var targetCount = _currentInstance.select({
            target: elInfo.id,
          }).length;
          return (
            elInfo.def.maxConnections > 0 &&
            targetCount >= elInfo.def.maxConnections
          );
        },
        element: elInfo.el,
        elementId: elInfo.id,
        isSource,
        isTarget,
        addClass: function (clazz) {
          _currentInstance.addClass(elInfo.el, clazz);
        },
        removeClass: function (clazz) {
          _currentInstance.removeClass(elInfo.el, clazz);
        },
        onDrop: function (jpc) {
          var source = jpc.endpoints[0];
          source.anchor.locked = false;
        },
        isDropAllowed: function () {
          return proxyComponent.isDropAllowed.apply(proxyComponent, arguments);
        },
        isRedrop: function (jpc) {
          return (
            jpc.suspendedElement != null &&
            jpc.suspendedEndpoint != null &&
            jpc.suspendedEndpoint.element === elInfo.el
          );
        },
        getEndpoint: function (jpc) {
          var newEndpoint = elInfo.def.endpoint;
          if (newEndpoint == null || newEndpoint._jsPlumb == null) {
            var eps = _currentInstance.deriveEndpointAndAnchorSpec(
              jpc.getType().join(" "),
              true
            );
            var pp = eps.endpoints
              ? root.jsPlumb.extend(p, {
                  endpoint: elInfo.def.def.endpoint || eps.endpoints[1],
                })
              : p;
            if (eps.anchors)
              pp = root.jsPlumb.extend(pp, {
                anchor: elInfo.def.def.anchor || eps.anchors[1],
              });
            newEndpoint = _currentInstance.addEndpoint(elInfo.el, pp);
            newEndpoint._mtNew = true;
          }
          if (p.uniqueEndpoint) elInfo.def.endpoint = newEndpoint;
          newEndpoint.setDeleteOnEmpty(true);
          if (jpc.isDetachable()) newEndpoint.initDraggable();
          if (newEndpoint.anchor.positionFinder != null) {
            var dropPosition = _currentInstance.getUIPosition(
              arguments,
              _currentInstance.getZoom()
            );
            var elPosition = _currentInstance.getOffset(elInfo.el);
            var elSize = _currentInstance.getSize(elInfo.el);
            var ap =
              dropPosition == null
                ? [0, 0]
                : newEndpoint.anchor.positionFinder(
                    dropPosition,
                    elPosition,
                    elSize,
                    newEndpoint.anchor.constructorParams
                  );
            newEndpoint.anchor.x = ap[0];
            newEndpoint.anchor.y = ap[1];
          }
          return newEndpoint;
        },
        maybeCleanup: function (ep) {
          if (ep._mtNew && ep.connections.length === 0)
            _currentInstance.deleteObject({ endpoint: ep });
          else delete ep._mtNew;
        },
      });
      var dropEvent = root.jsPlumb.dragEvents.drop;
      dropOptions.scope =
        dropOptions.scope || p.scope || _currentInstance.Defaults.Scope;
      dropOptions[dropEvent] = _ju.wrap(dropOptions[dropEvent], _drop, true);
      dropOptions.rank = p.rank || 0;
      if (isTarget)
        dropOptions[root.jsPlumb.dragEvents.over] = function () {
          return true;
        };
      if (p.allowLoopback === false)
        dropOptions.canDrop = function (_drag) {
          var de = _drag.getDragElement()._jsPlumbRelatedElement;
          return de !== elInfo.el;
        };
      _currentInstance.initDroppable(elInfo.el, dropOptions, "internal");
      return _drop;
    };
    this.makeTarget = function (el$jscomp$0, params, referenceParams) {
      var p = root.jsPlumb.extend({ _jsPlumb: this }, referenceParams);
      root.jsPlumb.extend(p, params);
      var maxConnections = p.maxConnections || -1;
      var _doOne = function (el) {
        var elInfo = _info(el);
        var elid = elInfo.id;
        var dropOptions = root.jsPlumb.extend({}, p.dropOptions || {});
        var type = p.connectionType || "default";
        this.targetEndpointDefinitions[elid] =
          this.targetEndpointDefinitions[elid] || {};
        _ensureContainer(elid);
        if (elInfo.el._isJsPlumbGroup && dropOptions.rank == null)
          dropOptions.rank = -1;
        var _def = {
          def: root.jsPlumb.extend({}, p),
          uniqueEndpoint: p.uniqueEndpoint,
          maxConnections,
          enabled: true,
        };
        if (p.createEndpoint) {
          _def.uniqueEndpoint = true;
          _def.endpoint = _currentInstance.addEndpoint(el, _def.def);
          _def.endpoint.setDeleteOnEmpty(false);
        }
        elInfo.def = _def;
        this.targetEndpointDefinitions[elid][type] = _def;
        _makeElementDropHandler(
          elInfo,
          p,
          dropOptions,
          p.isSource === true,
          true
        );
        elInfo.el._katavorioDrop[
          elInfo.el._katavorioDrop.length - 1
        ].targetDef = _def;
      }.bind(this);
      var inputs =
        el$jscomp$0.length && el$jscomp$0.constructor !== String
          ? el$jscomp$0
          : [el$jscomp$0];
      var i = 0;
      for (var ii = inputs.length; i < ii; i++) _doOne(inputs[i]);
      return this;
    };
    this.unmakeTarget = function (el, doNotClearArrays) {
      var info = _info(el);
      _currentInstance.destroyDroppable(info.el, "internal");
      if (!doNotClearArrays) delete this.targetEndpointDefinitions[info.id];
      return this;
    };
    this.makeSource = function (el, params, referenceParams) {
      var p = root.jsPlumb.extend({ _jsPlumb: this }, referenceParams);
      root.jsPlumb.extend(p, params);
      var type = p.connectionType || "default";
      var aae = _currentInstance.deriveEndpointAndAnchorSpec(type);
      p.endpoint = p.endpoint || aae.endpoints[0];
      p.anchor = p.anchor || aae.anchors[0];
      var maxConnections = p.maxConnections || -1;
      var onMaxConnections = p.onMaxConnections;
      var _doOne = function (elInfo) {
        var elid = elInfo.id;
        var _del = this.getElement(elInfo.el);
        this.sourceEndpointDefinitions[elid] =
          this.sourceEndpointDefinitions[elid] || {};
        _ensureContainer(elid);
        var _def = {
          def: root.jsPlumb.extend({}, p),
          uniqueEndpoint: p.uniqueEndpoint,
          maxConnections,
          enabled: true,
        };
        if (p.createEndpoint) {
          _def.uniqueEndpoint = true;
          _def.endpoint = _currentInstance.addEndpoint(el, _def.def);
          _def.endpoint.setDeleteOnEmpty(false);
        }
        this.sourceEndpointDefinitions[elid][type] = _def;
        elInfo.def = _def;
        var stopEvent = root.jsPlumb.dragEvents.stop;
        var dragEvent = root.jsPlumb.dragEvents.drag;
        var dragOptions = root.jsPlumb.extend({}, p.dragOptions || {});
        var existingDrag = dragOptions.drag;
        var existingStop = dragOptions.stop;
        var ep = null;
        var endpointAddedButNoDragYet = false;
        dragOptions.scope = dragOptions.scope || p.scope;
        dragOptions[dragEvent] = _ju.wrap(dragOptions[dragEvent], function () {
          if (existingDrag) existingDrag.apply(this, arguments);
          endpointAddedButNoDragYet = false;
        });
        dragOptions[stopEvent] = _ju.wrap(
          dragOptions[stopEvent],
          function () {
            if (existingStop) existingStop.apply(this, arguments);
            this.currentlyDragging = false;
            if (ep._jsPlumb != null) {
              var anchorDef = p.anchor || this.Defaults.Anchor;
              var oldAnchor = ep.anchor;
              var oldConnection = ep.connections[0];
              var newAnchor = this.makeAnchor(anchorDef, elid, this);
              var _el = ep.element;
              if (newAnchor.positionFinder != null) {
                var elPosition = _currentInstance.getOffset(_el);
                var elSize = this.getSize(_el);
                var dropPosition = {
                  left: elPosition.left + oldAnchor.x * elSize[0],
                  top: elPosition.top + oldAnchor.y * elSize[1],
                };
                var ap = newAnchor.positionFinder(
                  dropPosition,
                  elPosition,
                  elSize,
                  newAnchor.constructorParams
                );
                newAnchor.x = ap[0];
                newAnchor.y = ap[1];
              }
              ep.setAnchor(newAnchor, true);
              ep.repaint();
              this.repaint(ep.elementId);
              if (oldConnection != null) this.repaint(oldConnection.targetId);
            }
          }.bind(this)
        );
        var mouseDownListener = function (e) {
          if (e.which === 3 || e.button === 2) return;
          elid = this.getId(this.getElement(elInfo.el));
          var def = this.sourceEndpointDefinitions[elid][type];
          if (!def.enabled) return;
          if (p.filter) {
            var r = _ju.isString(p.filter)
              ? selectorFilter(e, elInfo.el, p.filter, this, p.filterExclude)
              : p.filter(e, elInfo.el);
            if (r === false) return;
          }
          var sourceCount = this.select({ source: elid }).length;
          if (def.maxConnections >= 0 && sourceCount >= def.maxConnections) {
            if (onMaxConnections)
              onMaxConnections({ element: elInfo.el, maxConnections }, e);
            return false;
          }
          var elxy = root.jsPlumb.getPositionOnElement(e, _del, _zoom);
          var tempEndpointParams = {};
          root.jsPlumb.extend(tempEndpointParams, def.def);
          tempEndpointParams.isTemporarySource = true;
          tempEndpointParams.anchor = [elxy[0], elxy[1], 0, 0];
          tempEndpointParams.dragOptions = dragOptions;
          if (def.def.scope) tempEndpointParams.scope = def.def.scope;
          ep = this.addEndpoint(elid, tempEndpointParams);
          endpointAddedButNoDragYet = true;
          ep.setDeleteOnEmpty(true);
          if (def.uniqueEndpoint)
            if (!def.endpoint) {
              def.endpoint = ep;
              ep.setDeleteOnEmpty(false);
            } else ep.finalEndpoint = def.endpoint;
          var _delTempEndpoint = function () {
            _currentInstance.off(ep.canvas, "mouseup", _delTempEndpoint);
            _currentInstance.off(elInfo.el, "mouseup", _delTempEndpoint);
            if (endpointAddedButNoDragYet) {
              endpointAddedButNoDragYet = false;
              _currentInstance.deleteEndpoint(ep);
            }
          };
          _currentInstance.on(ep.canvas, "mouseup", _delTempEndpoint);
          _currentInstance.on(elInfo.el, "mouseup", _delTempEndpoint);
          var payload = {};
          if (def.def.extract)
            for (var att in def.def.extract) {
              var v = (e.srcElement || e.target).getAttribute(att);
              if (v) payload[def.def.extract[att]] = v;
            }
          _currentInstance.trigger(ep.canvas, "mousedown", e, payload);
          _ju.consume(e);
        }.bind(this);
        this.on(elInfo.el, "mousedown", mouseDownListener);
        _def.trigger = mouseDownListener;
        if (p.filter && (_ju.isString(p.filter) || _ju.isFunction(p.filter)))
          _currentInstance.setDragFilter(elInfo.el, p.filter);
        var dropOptions = root.jsPlumb.extend({}, p.dropOptions || {});
        _makeElementDropHandler(
          elInfo,
          p,
          dropOptions,
          true,
          p.isTarget === true
        );
      }.bind(this);
      var inputs = el.length && el.constructor !== String ? el : [el];
      var i = 0;
      for (var ii = inputs.length; i < ii; i++) _doOne(_info(inputs[i]));
      return this;
    };
    this.unmakeSource = function (el, connectionType, doNotClearArrays) {
      var info = _info(el);
      _currentInstance.destroyDroppable(info.el, "internal");
      var eldefs = this.sourceEndpointDefinitions[info.id];
      if (eldefs)
        for (var def in eldefs)
          if (connectionType == null || connectionType === def) {
            var mouseDownListener = eldefs[def].trigger;
            if (mouseDownListener)
              _currentInstance.off(info.el, "mousedown", mouseDownListener);
            if (!doNotClearArrays)
              delete this.sourceEndpointDefinitions[info.id][def];
          }
      return this;
    };
    this.unmakeEverySource = function () {
      for (var i in this.sourceEndpointDefinitions)
        _currentInstance.unmakeSource(i, null, true);
      this.sourceEndpointDefinitions = {};
      return this;
    };
    var _getScope = function (el, types, connectionType) {
      types = _ju.isArray(types) ? types : [types];
      var id = _getId(el);
      connectionType = connectionType || "default";
      for (var i = 0; i < types.length; i++) {
        var eldefs = this[types[i]][id];
        if (eldefs && eldefs[connectionType])
          return eldefs[connectionType].def.scope || this.Defaults.Scope;
      }
    }.bind(this);
    var _setScope = function (el, scope, types, connectionType) {
      types = _ju.isArray(types) ? types : [types];
      var id = _getId(el);
      connectionType = connectionType || "default";
      for (var i = 0; i < types.length; i++) {
        var eldefs = this[types[i]][id];
        if (eldefs && eldefs[connectionType])
          eldefs[connectionType].def.scope = scope;
      }
    }.bind(this);
    this.getScope = function (el, scope) {
      return _getScope(el, [
        "sourceEndpointDefinitions",
        "targetEndpointDefinitions",
      ]);
    };
    this.getSourceScope = function (el) {
      return _getScope(el, "sourceEndpointDefinitions");
    };
    this.getTargetScope = function (el) {
      return _getScope(el, "targetEndpointDefinitions");
    };
    this.setScope = function (el, scope, connectionType) {
      this.setSourceScope(el, scope, connectionType);
      this.setTargetScope(el, scope, connectionType);
    };
    this.setSourceScope = function (el, scope, connectionType) {
      _setScope(el, scope, "sourceEndpointDefinitions", connectionType);
      this.setDragScope(el, scope);
    };
    this.setTargetScope = function (el, scope, connectionType) {
      _setScope(el, scope, "targetEndpointDefinitions", connectionType);
      this.setDropScope(el, scope);
    };
    this.unmakeEveryTarget = function () {
      for (var i in this.targetEndpointDefinitions)
        _currentInstance.unmakeTarget(i, true);
      this.targetEndpointDefinitions = {};
      return this;
    };
    var _setEnabled = function (type, el, state, toggle, connectionType) {
      var a =
        type === "source"
          ? this.sourceEndpointDefinitions
          : this.targetEndpointDefinitions;
      connectionType = connectionType || "default";
      if (el.length && !_ju.isString(el)) {
        var originalState = [];
        var i = 0;
        for (var ii = el.length; i < ii; i++) {
          var info = _info(el[i]);
          if (a[info.id] && a[info.id][connectionType]) {
            originalState[i] = a[info.id][connectionType].enabled;
            var newState = toggle ? !originalState[i] : state;
            a[info.id][connectionType].enabled = newState;
            _currentInstance[newState ? "removeClass" : "addClass"](
              info.el,
              "jtk-" + type + "-disabled"
            );
          }
        }
      } else {
        info = _info(el);
        var id = info.id;
        if (a[id] && a[id][connectionType]) {
          originalState = a[id][connectionType].enabled;
          newState = toggle ? !originalState : state;
          a[id][connectionType].enabled = newState;
          _currentInstance[newState ? "removeClass" : "addClass"](
            info.el,
            "jtk-" + type + "-disabled"
          );
        }
      }
      return originalState;
    }.bind(this);
    var _first = function (el, fn) {
      if (el != null)
        if (_ju.isString(el) || !el.length) return fn.apply(this, [el]);
        else if (el.length) return fn.apply(this, [el[0]]);
    }.bind(this);
    this.toggleSourceEnabled = function (el, connectionType) {
      _setEnabled("source", el, null, true, connectionType);
      return this.isSourceEnabled(el, connectionType);
    };
    this.setSourceEnabled = function (el, state, connectionType) {
      return _setEnabled("source", el, state, null, connectionType);
    };
    this.isSource = function (el, connectionType) {
      connectionType = connectionType || "default";
      return _first(
        el,
        function (_el) {
          var eldefs = this.sourceEndpointDefinitions[_info(_el).id];
          return eldefs != null && eldefs[connectionType] != null;
        }.bind(this)
      );
    };
    this.isSourceEnabled = function (el, connectionType) {
      connectionType = connectionType || "default";
      return _first(
        el,
        function (_el) {
          var sep = this.sourceEndpointDefinitions[_info(_el).id];
          return (
            sep && sep[connectionType] && sep[connectionType].enabled === true
          );
        }.bind(this)
      );
    };
    this.toggleTargetEnabled = function (el, connectionType) {
      _setEnabled("target", el, null, true, connectionType);
      return this.isTargetEnabled(el, connectionType);
    };
    this.isTarget = function (el, connectionType) {
      connectionType = connectionType || "default";
      return _first(
        el,
        function (_el) {
          var eldefs = this.targetEndpointDefinitions[_info(_el).id];
          return eldefs != null && eldefs[connectionType] != null;
        }.bind(this)
      );
    };
    this.isTargetEnabled = function (el, connectionType) {
      connectionType = connectionType || "default";
      return _first(
        el,
        function (_el) {
          var tep = this.targetEndpointDefinitions[_info(_el).id];
          return (
            tep && tep[connectionType] && tep[connectionType].enabled === true
          );
        }.bind(this)
      );
    };
    this.setTargetEnabled = function (el, state, connectionType) {
      return _setEnabled("target", el, state, null, connectionType);
    };
    this.ready = function (fn) {
      _currentInstance.bind("ready", fn);
    };
    var _elEach = function (el, fn) {
      if (typeof el === "object" && el.length) {
        var i = 0;
        for (var ii = el.length; i < ii; i++) fn(el[i]);
      } else fn(el);
      return _currentInstance;
    };
    this.repaint = function (el, ui, timestamp) {
      return _elEach(el, function (_el) {
        _draw(_el, ui, timestamp);
      });
    };
    this.revalidate = function (el, timestamp, isIdAlready) {
      var elId = isIdAlready ? el : _currentInstance.getId(el);
      _currentInstance.updateOffset({ elId, recalc: true, timestamp });
      var dm = _currentInstance.getDragManager();
      if (dm) dm.updateOffsets(elId);
      return _draw(el, null, timestamp);
    };
    this.repaintEverything = function () {
      var timestamp = jsPlumbUtil.uuid();
      for (var elId in endpointsByElement)
        _currentInstance.updateOffset({ elId, recalc: true, timestamp });
      for (elId in endpointsByElement) _draw(elId, null, timestamp);
      return this;
    };
    this.removeAllEndpoints = function (el, recurse, affectedElements) {
      affectedElements = affectedElements || [];
      var _one = function (_el) {
        var info = _info(_el);
        var ebe = endpointsByElement[info.id];
        var i;
        var ii;
        if (ebe) {
          affectedElements.push(info);
          for (i = 0, ii = ebe.length; i < ii; i++)
            _currentInstance.deleteEndpoint(ebe[i], false);
        }
        delete endpointsByElement[info.id];
        if (recurse)
          if (info.el && info.el.nodeType !== 3 && info.el.nodeType !== 8)
            for (i = 0, ii = info.el.childNodes.length; i < ii; i++)
              _one(info.el.childNodes[i]);
      };
      _one(el);
      return this;
    };
    var _doRemove = function (info, affectedElements) {
      _currentInstance.removeAllEndpoints(info.id, true, affectedElements);
      var dm = _currentInstance.getDragManager();
      var _one = function (_info) {
        if (dm) dm.elementRemoved(_info.id);
        _currentInstance.router.elementRemoved(_info.id);
        if (_currentInstance.isSource(_info.el))
          _currentInstance.unmakeSource(_info.el);
        if (_currentInstance.isTarget(_info.el))
          _currentInstance.unmakeTarget(_info.el);
        _currentInstance.destroyDraggable(_info.el);
        _currentInstance.destroyDroppable(_info.el);
        delete _currentInstance.floatingConnections[_info.id];
        delete managedElements[_info.id];
        delete offsets[_info.id];
        if (_info.el) {
          _currentInstance.removeElement(_info.el);
          _info.el._jsPlumb = null;
        }
      };
      for (var ae = 1; ae < affectedElements.length; ae++)
        _one(affectedElements[ae]);
      _one(info);
    };
    this.remove = function (el, doNotRepaint) {
      var info = _info(el);
      var affectedElements = [];
      if (info.text && info.el.parentNode)
        info.el.parentNode.removeChild(info.el);
      else if (info.id)
        _currentInstance.batch(function () {
          _doRemove(info, affectedElements);
        }, doNotRepaint === true);
      return _currentInstance;
    };
    this.empty = function (el$jscomp$0, doNotRepaint) {
      var affectedElements = [];
      var _one = function (el, dontRemoveFocus) {
        var info = _info(el);
        if (info.text) info.el.parentNode.removeChild(info.el);
        else if (info.el) {
          for (; info.el.childNodes.length > 0; ) _one(info.el.childNodes[0]);
          if (!dontRemoveFocus) _doRemove(info, affectedElements);
        }
      };
      _currentInstance.batch(function () {
        _one(el$jscomp$0, true);
      }, doNotRepaint === false);
      return _currentInstance;
    };
    this.reset = function (doNotUnbindInstanceEventListeners) {
      _currentInstance.silently(
        function () {
          _hoverSuspended = false;
          _currentInstance.removeAllGroups();
          _currentInstance.removeGroupManager();
          _currentInstance.deleteEveryEndpoint();
          if (!doNotUnbindInstanceEventListeners) _currentInstance.unbind();
          this.targetEndpointDefinitions = {};
          this.sourceEndpointDefinitions = {};
          connections.length = 0;
          if (this.doReset) this.doReset();
        }.bind(this)
      );
    };
    this.destroy = function () {
      this.reset();
      _container = null;
      _containerDelegations = null;
    };
    var _clearObject = function (obj) {
      if (obj.canvas && obj.canvas.parentNode)
        obj.canvas.parentNode.removeChild(obj.canvas);
      obj.cleanup();
      obj.destroy();
    };
    this.clear = function () {
      _currentInstance.select().each(_clearObject);
      _currentInstance.selectEndpoints().each(_clearObject);
      endpointsByElement = {};
      endpointsByUUID = {};
    };
    this.setDefaultScope = function (scope) {
      DEFAULT_SCOPE = scope;
      return _currentInstance;
    };
    this.deriveEndpointAndAnchorSpec = function (type, dontPrependDefault) {
      var bits = ((dontPrependDefault ? "" : "default ") + type).split(/[\s]/);
      var eps = null;
      var ep = null;
      var a = null;
      var as = null;
      for (var i = 0; i < bits.length; i++) {
        var _t = _currentInstance.getType(bits[i], "connection");
        if (_t) {
          if (_t.endpoints) eps = _t.endpoints;
          if (_t.endpoint) ep = _t.endpoint;
          if (_t.anchors) as = _t.anchors;
          if (_t.anchor) a = _t.anchor;
        }
      }
      return { endpoints: eps ? eps : [ep, ep], anchors: as ? as : [a, a] };
    };
    this.setId = function (el, newId, doNotSetAttribute) {
      if (_ju.isString(el)) var id = el;
      else {
        el = this.getElement(el);
        id = this.getId(el);
      }
      var sConns = this.getConnections({ source: id, scope: "*" }, true);
      var tConns = this.getConnections({ target: id, scope: "*" }, true);
      newId = "" + newId;
      if (!doNotSetAttribute) {
        el = this.getElement(id);
        this.setAttribute(el, "id", newId);
      } else el = this.getElement(newId);
      endpointsByElement[newId] = endpointsByElement[id] || [];
      var i$jscomp$0 = 0;
      for (
        var ii$jscomp$0 = endpointsByElement[newId].length;
        i$jscomp$0 < ii$jscomp$0;
        i$jscomp$0++
      ) {
        endpointsByElement[newId][i$jscomp$0].setElementId(newId);
        endpointsByElement[newId][i$jscomp$0].setReferenceElement(el);
      }
      delete endpointsByElement[id];
      this.sourceEndpointDefinitions[newId] =
        this.sourceEndpointDefinitions[id];
      delete this.sourceEndpointDefinitions[id];
      this.targetEndpointDefinitions[newId] =
        this.targetEndpointDefinitions[id];
      delete this.targetEndpointDefinitions[id];
      this.router.changeId(id, newId);
      var dm = this.getDragManager();
      if (dm) dm.changeId(id, newId);
      managedElements[newId] = managedElements[id];
      delete managedElements[id];
      var _conns = function (list, epIdx, type) {
        var i = 0;
        for (var ii = list.length; i < ii; i++) {
          list[i].endpoints[epIdx].setElementId(newId);
          list[i].endpoints[epIdx].setReferenceElement(el);
          list[i][type + "Id"] = newId;
          list[i][type] = el;
        }
      };
      _conns(sConns, 0, "source");
      _conns(tConns, 1, "target");
      this.repaint(newId);
    };
    this.setDebugLog = function (debugLog) {
      log = debugLog;
    };
    this.setSuspendDrawing = function (val, repaintAfterwards) {
      var curVal = _suspendDrawing;
      _suspendDrawing = val;
      if (val) _suspendedAt = new Date().getTime();
      else _suspendedAt = null;
      if (repaintAfterwards) this.repaintEverything();
      return curVal;
    };
    this.isSuspendDrawing = function () {
      return _suspendDrawing;
    };
    this.getSuspendedAt = function () {
      return _suspendedAt;
    };
    this.batch = function (fn, doNotRepaintAfterwards) {
      var _wasSuspended = this.isSuspendDrawing();
      if (!_wasSuspended) this.setSuspendDrawing(true);
      try {
        fn();
      } catch (e) {
        _ju.log("Function run while suspended failed", e);
      }
      if (!_wasSuspended)
        this.setSuspendDrawing(false, !doNotRepaintAfterwards);
    };
    this.doWhileSuspended = this.batch;
    this.getCachedData = _getCachedData;
    this.show = function (el, changeEndpoints) {
      _setVisible(el, "block", changeEndpoints);
      return _currentInstance;
    };
    this.toggleVisible = _toggleVisible;
    this.addListener = this.bind;
    var floatingConnections = [];
    this.registerFloatingConnection = function (info, conn, ep) {
      floatingConnections[info.id] = conn;
      _ju.addToList(endpointsByElement, info.id, ep);
    };
    this.getFloatingConnectionFor = function (id) {
      return floatingConnections[id];
    };
    this.listManager = new root.jsPlumbListManager(
      this,
      this.Defaults.ListStyle
    );
  });
  _ju.extend(root.jsPlumbInstance, _ju.EventGenerator, {
    setAttribute: function (el, a, v) {
      this.setAttribute(el, a, v);
    },
    getAttribute: function (el, a) {
      return this.getAttribute(root.jsPlumb.getElement(el), a);
    },
    convertToFullOverlaySpec: function (spec) {
      if (_ju.isString(spec)) spec = [spec, {}];
      spec[1].id = spec[1].id || _ju.uuid();
      return spec;
    },
    registerConnectionType: function (id, type) {
      this._connectionTypes[id] = root.jsPlumb.extend({}, type);
      if (type.overlays) {
        var to = {};
        for (var i = 0; i < type.overlays.length; i++) {
          var fo = this.convertToFullOverlaySpec(type.overlays[i]);
          to[fo[1].id] = fo;
        }
        this._connectionTypes[id].overlays = to;
      }
    },
    registerConnectionTypes: function (types) {
      for (var i in types) this.registerConnectionType(i, types[i]);
    },
    registerEndpointType: function (id, type) {
      this._endpointTypes[id] = root.jsPlumb.extend({}, type);
      if (type.overlays) {
        var to = {};
        for (var i = 0; i < type.overlays.length; i++) {
          var fo = this.convertToFullOverlaySpec(type.overlays[i]);
          to[fo[1].id] = fo;
        }
        this._endpointTypes[id].overlays = to;
      }
    },
    registerEndpointTypes: function (types) {
      for (var i in types) this.registerEndpointType(i, types[i]);
    },
    getType: function (id, typeDescriptor) {
      return typeDescriptor === "connection"
        ? this._connectionTypes[id]
        : this._endpointTypes[id];
    },
    setIdChanged: function (oldId, newId) {
      this.setId(oldId, newId, true);
    },
    setParent: function (el, newParent) {
      var _dom = this.getElement(el);
      var _id = this.getId(_dom);
      var _pdom = this.getElement(newParent);
      var _pid = this.getId(_pdom);
      var dm = this.getDragManager();
      _dom.parentNode.removeChild(_dom);
      _pdom.appendChild(_dom);
      if (dm) dm.setParent(_dom, _id, _pdom, _pid);
    },
    extend: function (o1, o2, names) {
      var i;
      if (names) for (i = 0; i < names.length; i++) o1[names[i]] = o2[names[i]];
      else for (i in o2) o1[i] = o2[i];
      return o1;
    },
    floatingConnections: {},
    getFloatingAnchorIndex: function (jpc) {
      return jpc.endpoints[0].isFloating()
        ? 0
        : jpc.endpoints[1].isFloating()
        ? 1
        : -1;
    },
    proxyConnection: function (
      connection,
      index,
      proxyEl,
      proxyElId,
      endpointGenerator,
      anchorGenerator
    ) {
      var originalElementId = connection.endpoints[index].elementId;
      var originalEndpoint = connection.endpoints[index];
      connection.proxies = connection.proxies || [];
      if (connection.proxies[index]) var proxyEp = connection.proxies[index].ep;
      else
        proxyEp = this.addEndpoint(proxyEl, {
          endpoint: endpointGenerator(connection, index),
          anchor: anchorGenerator(connection, index),
          parameters: { isProxyEndpoint: true },
        });
      proxyEp.setDeleteOnEmpty(true);
      connection.proxies[index] = { ep: proxyEp, originalEp: originalEndpoint };
      if (index === 0)
        this.router.sourceOrTargetChanged(
          originalElementId,
          proxyElId,
          connection,
          proxyEl,
          0
        );
      else
        this.router.sourceOrTargetChanged(
          originalElementId,
          proxyElId,
          connection,
          proxyEl,
          1
        );
      originalEndpoint.detachFromConnection(connection, null, true);
      proxyEp.connections = [connection];
      connection.endpoints[index] = proxyEp;
      originalEndpoint.setVisible(false);
      connection.setVisible(true);
      this.revalidate(proxyEl);
    },
    unproxyConnection: function (connection, index, proxyElId) {
      if (
        connection._jsPlumb == null ||
        connection.proxies == null ||
        connection.proxies[index] == null
      )
        return;
      var originalElement = connection.proxies[index].originalEp.element;
      var originalElementId = connection.proxies[index].originalEp.elementId;
      connection.endpoints[index] = connection.proxies[index].originalEp;
      if (index === 0)
        this.router.sourceOrTargetChanged(
          proxyElId,
          originalElementId,
          connection,
          originalElement,
          0
        );
      else
        this.router.sourceOrTargetChanged(
          proxyElId,
          originalElementId,
          connection,
          originalElement,
          1
        );
      connection.proxies[index].ep.detachFromConnection(connection, null);
      connection.proxies[index].originalEp.addConnection(connection);
      if (connection.isVisible())
        connection.proxies[index].originalEp.setVisible(true);
      delete connection.proxies[index];
    },
  });
  var jsPlumb = new jsPlumbInstance$jscomp$0();
  root.jsPlumb = jsPlumb;
  jsPlumb.getInstance = function (_defaults, overrideFns) {
    var j = new jsPlumbInstance$jscomp$0(_defaults);
    if (overrideFns) for (var ovf in overrideFns) j[ovf] = overrideFns[ovf];
    j.init();
    return j;
  };
  jsPlumb.each = function (spec, fn) {
    if (spec == null) return;
    if (typeof spec === "string") fn(jsPlumb.getElement(spec));
    else if (spec.length != null)
      for (var i = 0; i < spec.length; i++) fn(jsPlumb.getElement(spec[i]));
    else fn(spec);
  };
  if (typeof exports !== "undefined") exports.jsPlumb = jsPlumb;
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var _internalLabelOverlayId = "__label";
  var _makeLabelOverlay = function (component, params) {
    var _params = {
      cssClass: params.cssClass,
      labelStyle: component.labelStyle,
      id: _internalLabelOverlayId,
      component,
      _jsPlumb: component._jsPlumb.instance,
    };
    var mergedParams = _jp.extend(_params, params);
    return new _jp.Overlays[component._jsPlumb.instance.getRenderMode()].Label(
      mergedParams
    );
  };
  var _processOverlay = function (component, o) {
    var _newOverlay = null;
    if (_ju.isArray(o)) {
      var type = o[0];
      var p = _jp.extend(
        { component, _jsPlumb: component._jsPlumb.instance },
        o[1]
      );
      if (o.length === 3) _jp.extend(p, o[2]);
      _newOverlay = new _jp.Overlays[
        component._jsPlumb.instance.getRenderMode()
      ][type](p);
    } else if (o.constructor === String)
      _newOverlay = new _jp.Overlays[
        component._jsPlumb.instance.getRenderMode()
      ][o]({ component, _jsPlumb: component._jsPlumb.instance });
    else _newOverlay = o;
    _newOverlay.id = _newOverlay.id || _ju.uuid();
    component.cacheTypeItem("overlay", _newOverlay, _newOverlay.id);
    component._jsPlumb.overlays[_newOverlay.id] = _newOverlay;
    return _newOverlay;
  };
  _jp.OverlayCapableJsPlumbUIComponent = function (params) {
    root.jsPlumbUIComponent.apply(this, arguments);
    this._jsPlumb.overlays = {};
    this._jsPlumb.overlayPositions = {};
    if (params.label)
      this.getDefaultType().overlays[_internalLabelOverlayId] = [
        "Label",
        {
          label: params.label,
          location: params.labelLocation || this.defaultLabelLocation || 0.5,
          labelStyle:
            params.labelStyle || this._jsPlumb.instance.Defaults.LabelStyle,
          id: _internalLabelOverlayId,
        },
      ];
    this.setListenerComponent = function (c) {
      if (this._jsPlumb)
        for (var i in this._jsPlumb.overlays)
          this._jsPlumb.overlays[i].setListenerComponent(c);
    };
  };
  _jp.OverlayCapableJsPlumbUIComponent.applyType = function (component, t) {
    if (t.overlays) {
      var keep = {};
      for (var i in t.overlays) {
        var existing = component._jsPlumb.overlays[t.overlays[i][1].id];
        if (existing) {
          existing.updateFrom(t.overlays[i][1]);
          keep[t.overlays[i][1].id] = true;
          existing.reattach(component._jsPlumb.instance, component);
        } else {
          var c = component.getCachedTypeItem("overlay", t.overlays[i][1].id);
          if (c != null) {
            c.reattach(component._jsPlumb.instance, component);
            c.setVisible(true);
            c.updateFrom(t.overlays[i][1]);
            component._jsPlumb.overlays[c.id] = c;
          } else c = component.addOverlay(t.overlays[i], true);
          keep[c.id] = true;
        }
      }
      for (i in component._jsPlumb.overlays)
        if (keep[component._jsPlumb.overlays[i].id] == null)
          component.removeOverlay(component._jsPlumb.overlays[i].id, true);
    }
  };
  _ju.extend(_jp.OverlayCapableJsPlumbUIComponent, root.jsPlumbUIComponent, {
    setHover: function (hover, ignoreAttachedElements) {
      if (this._jsPlumb && !this._jsPlumb.instance.isConnectionBeingDragged())
        for (var i in this._jsPlumb.overlays)
          this._jsPlumb.overlays[i][hover ? "addClass" : "removeClass"](
            this._jsPlumb.instance.hoverClass
          );
    },
    addOverlay: function (overlay, doNotRepaint) {
      var o = _processOverlay(this, overlay);
      if (this.getData && o.type === "Label" && _ju.isArray(overlay)) {
        var d = this.getData();
        var p = overlay[1];
        if (d) {
          var locationAttribute = p.labelLocationAttribute || "labelLocation";
          var loc = d ? d[locationAttribute] : null;
          if (loc) o.loc = loc;
        }
      }
      if (!doNotRepaint) this.repaint();
      return o;
    },
    getOverlay: function (id) {
      return this._jsPlumb.overlays[id];
    },
    getOverlays: function () {
      return this._jsPlumb.overlays;
    },
    hideOverlay: function (id) {
      var o = this.getOverlay(id);
      if (o) o.hide();
    },
    hideOverlays: function () {
      for (var i in this._jsPlumb.overlays) this._jsPlumb.overlays[i].hide();
    },
    showOverlay: function (id) {
      var o = this.getOverlay(id);
      if (o) o.show();
    },
    showOverlays: function () {
      for (var i in this._jsPlumb.overlays) this._jsPlumb.overlays[i].show();
    },
    removeAllOverlays: function (doNotRepaint) {
      for (var i in this._jsPlumb.overlays)
        if (this._jsPlumb.overlays[i].cleanup)
          this._jsPlumb.overlays[i].cleanup();
      this._jsPlumb.overlays = {};
      this._jsPlumb.overlayPositions = null;
      this._jsPlumb.overlayPlacements = {};
      if (!doNotRepaint) this.repaint();
    },
    removeOverlay: function (overlayId, dontCleanup) {
      var o = this._jsPlumb.overlays[overlayId];
      if (o) {
        o.setVisible(false);
        if (!dontCleanup && o.cleanup) o.cleanup();
        delete this._jsPlumb.overlays[overlayId];
        if (this._jsPlumb.overlayPositions)
          delete this._jsPlumb.overlayPositions[overlayId];
        if (this._jsPlumb.overlayPlacements)
          delete this._jsPlumb.overlayPlacements[overlayId];
      }
    },
    removeOverlays: function () {
      var i = 0;
      for (var j = arguments.length; i < j; i++)
        this.removeOverlay(arguments[i]);
    },
    moveParent: function (newParent) {
      if (this.bgCanvas) {
        this.bgCanvas.parentNode.removeChild(this.bgCanvas);
        newParent.appendChild(this.bgCanvas);
      }
      if (this.canvas && this.canvas.parentNode) {
        this.canvas.parentNode.removeChild(this.canvas);
        newParent.appendChild(this.canvas);
        for (var i in this._jsPlumb.overlays)
          if (this._jsPlumb.overlays[i].isAppendedAtTopLevel) {
            var el = this._jsPlumb.overlays[i].getElement();
            el.parentNode.removeChild(el);
            newParent.appendChild(el);
          }
      }
    },
    getLabel: function () {
      var lo = this.getOverlay(_internalLabelOverlayId);
      return lo != null ? lo.getLabel() : null;
    },
    getLabelOverlay: function () {
      return this.getOverlay(_internalLabelOverlayId);
    },
    setLabel: function (l) {
      var lo = this.getOverlay(_internalLabelOverlayId);
      if (!lo) {
        var params =
          l.constructor === String || l.constructor === Function
            ? { label: l }
            : l;
        lo = _makeLabelOverlay(this, params);
        this._jsPlumb.overlays[_internalLabelOverlayId] = lo;
      } else if (l.constructor === String || l.constructor === Function)
        lo.setLabel(l);
      else {
        if (l.label) lo.setLabel(l.label);
        if (l.location) lo.setLocation(l.location);
      }
      if (!this._jsPlumb.instance.isSuspendDrawing()) this.repaint();
    },
    cleanup: function (force) {
      for (var i in this._jsPlumb.overlays) {
        this._jsPlumb.overlays[i].cleanup(force);
        this._jsPlumb.overlays[i].destroy(force);
      }
      if (force) {
        this._jsPlumb.overlays = {};
        this._jsPlumb.overlayPositions = null;
      }
    },
    setVisible: function (v) {
      this[v ? "showOverlays" : "hideOverlays"]();
    },
    setAbsoluteOverlayPosition: function (overlay, xy) {
      this._jsPlumb.overlayPositions[overlay.id] = xy;
    },
    getAbsoluteOverlayPosition: function (overlay) {
      return this._jsPlumb.overlayPositions
        ? this._jsPlumb.overlayPositions[overlay.id]
        : null;
    },
    _clazzManip: function (action, clazz, dontUpdateOverlays) {
      if (!dontUpdateOverlays)
        for (var i in this._jsPlumb.overlays)
          this._jsPlumb.overlays[i][action + "Class"](clazz);
    },
    addClass: function (clazz, dontUpdateOverlays) {
      this._clazzManip("add", clazz, dontUpdateOverlays);
    },
    removeClass: function (clazz, dontUpdateOverlays) {
      this._clazzManip("remove", clazz, dontUpdateOverlays);
    },
  });
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var _makeConnectionDragHandler = function (endpoint, placeholder, _jsPlumb) {
    var stopped = false;
    return {
      drag: function () {
        if (stopped) {
          stopped = false;
          return true;
        }
        if (placeholder.element) {
          var _ui = _jsPlumb.getUIPosition(arguments, _jsPlumb.getZoom());
          if (_ui != null) _jsPlumb.setPosition(placeholder.element, _ui);
          _jsPlumb.repaint(placeholder.element, _ui);
          endpoint.paint({
            anchorPoint: endpoint.anchor.getCurrentLocation({
              element: endpoint,
            }),
          });
        }
      },
      stopDrag: function () {
        stopped = true;
      },
    };
  };
  var _makeDraggablePlaceholder = function (placeholder, _jsPlumb, ipco, ips) {
    var n = _jsPlumb.createElement("div", { position: "absolute" });
    _jsPlumb.appendElement(n);
    var id = _jsPlumb.getId(n);
    _jsPlumb.setPosition(n, ipco);
    n.style.width = ips[0] + "px";
    n.style.height = ips[1] + "px";
    _jsPlumb.manage(id, n, true);
    placeholder.id = id;
    placeholder.element = n;
  };
  var _makeFloatingEndpoint = function (
    paintStyle,
    referenceAnchor,
    endpoint,
    referenceCanvas,
    sourceElement,
    _jsPlumb,
    _newEndpoint,
    scope
  ) {
    var floatingAnchor = new _jp.FloatingAnchor({
      reference: referenceAnchor,
      referenceCanvas,
      jsPlumbInstance: _jsPlumb,
    });
    return _newEndpoint({
      paintStyle,
      endpoint,
      anchor: floatingAnchor,
      source: sourceElement,
      scope,
    });
  };
  var typeParameters = [
    "connectorStyle",
    "connectorHoverStyle",
    "connectorOverlays",
    "connector",
    "connectionType",
    "connectorClass",
    "connectorHoverClass",
  ];
  var findConnectionToUseForDynamicAnchor = function (
    ep,
    elementWithPrecedence
  ) {
    var idx = 0;
    if (elementWithPrecedence != null)
      for (var i = 0; i < ep.connections.length; i++)
        if (
          ep.connections[i].sourceId === elementWithPrecedence ||
          ep.connections[i].targetId === elementWithPrecedence
        ) {
          idx = i;
          break;
        }
    return ep.connections[idx];
  };
  _jp.Endpoint = function (params$jscomp$0) {
    var _jsPlumb = params$jscomp$0._jsPlumb;
    var _newConnection = params$jscomp$0.newConnection;
    var _newEndpoint = params$jscomp$0.newEndpoint;
    this.idPrefix = "_jsplumb_e_";
    this.defaultLabelLocation = [0.5, 0.5];
    this.defaultOverlayKeys = ["Overlays", "EndpointOverlays"];
    _jp.OverlayCapableJsPlumbUIComponent.apply(this, arguments);
    this.appendToDefaultType({
      connectionType: params$jscomp$0.connectionType,
      maxConnections:
        params$jscomp$0.maxConnections == null
          ? this._jsPlumb.instance.Defaults.MaxConnections
          : params$jscomp$0.maxConnections,
      paintStyle:
        params$jscomp$0.endpointStyle ||
        params$jscomp$0.paintStyle ||
        params$jscomp$0.style ||
        this._jsPlumb.instance.Defaults.EndpointStyle ||
        _jp.Defaults.EndpointStyle,
      hoverPaintStyle:
        params$jscomp$0.endpointHoverStyle ||
        params$jscomp$0.hoverPaintStyle ||
        this._jsPlumb.instance.Defaults.EndpointHoverStyle ||
        _jp.Defaults.EndpointHoverStyle,
      connectorStyle: params$jscomp$0.connectorStyle,
      connectorHoverStyle: params$jscomp$0.connectorHoverStyle,
      connectorClass: params$jscomp$0.connectorClass,
      connectorHoverClass: params$jscomp$0.connectorHoverClass,
      connectorOverlays: params$jscomp$0.connectorOverlays,
      connector: params$jscomp$0.connector,
      connectorTooltip: params$jscomp$0.connectorTooltip,
    });
    this._jsPlumb.enabled = !(params$jscomp$0.enabled === false);
    this._jsPlumb.visible = true;
    this.element = _jp.getElement(params$jscomp$0.source);
    this._jsPlumb.uuid = params$jscomp$0.uuid;
    this._jsPlumb.floatingEndpoint = null;
    var inPlaceCopy = null;
    if (this._jsPlumb.uuid)
      params$jscomp$0.endpointsByUUID[this._jsPlumb.uuid] = this;
    this.elementId = params$jscomp$0.elementId;
    this.dragProxy = params$jscomp$0.dragProxy;
    this._jsPlumb.connectionCost = params$jscomp$0.connectionCost;
    this._jsPlumb.connectionsDirected = params$jscomp$0.connectionsDirected;
    this._jsPlumb.currentAnchorClass = "";
    this._jsPlumb.events = {};
    var deleteOnEmpty = params$jscomp$0.deleteOnEmpty === true;
    this.setDeleteOnEmpty = function (d) {
      deleteOnEmpty = d;
    };
    var _updateAnchorClass = function () {
      var oldAnchorClass =
        _jsPlumb.endpointAnchorClassPrefix +
        "-" +
        this._jsPlumb.currentAnchorClass;
      this._jsPlumb.currentAnchorClass = this.anchor.getCssClass();
      var anchorClass =
        _jsPlumb.endpointAnchorClassPrefix +
        (this._jsPlumb.currentAnchorClass
          ? "-" + this._jsPlumb.currentAnchorClass
          : "");
      this.removeClass(oldAnchorClass);
      this.addClass(anchorClass);
      _jp.updateClasses(this.element, anchorClass, oldAnchorClass);
    }.bind(this);
    this.prepareAnchor = function (anchorParams) {
      var a = this._jsPlumb.instance.makeAnchor(
        anchorParams,
        this.elementId,
        _jsPlumb
      );
      a.bind(
        "anchorChanged",
        function (currentAnchor) {
          this.fire("anchorChanged", { endpoint: this, anchor: currentAnchor });
          _updateAnchorClass();
        }.bind(this)
      );
      return a;
    };
    this.setPreparedAnchor = function (anchor, doNotRepaint) {
      this._jsPlumb.instance.continuousAnchorFactory.clear(this.elementId);
      this.anchor = anchor;
      _updateAnchorClass();
      if (!doNotRepaint) this._jsPlumb.instance.repaint(this.elementId);
      return this;
    };
    this.setAnchor = function (anchorParams, doNotRepaint) {
      var a = this.prepareAnchor(anchorParams);
      this.setPreparedAnchor(a, doNotRepaint);
      return this;
    };
    var internalHover = function (state) {
      if (this.connections.length > 0)
        for (var i = 0; i < this.connections.length; i++)
          this.connections[i].setHover(state, false);
      else this.setHover(state);
    }.bind(this);
    this.bind("mouseover", function () {
      internalHover(true);
    });
    this.bind("mouseout", function () {
      internalHover(false);
    });
    if (!params$jscomp$0._transient)
      this._jsPlumb.instance.router.addEndpoint(this, this.elementId);
    this.prepareEndpoint = function (ep, typeId) {
      var _e = function (t, p) {
        var rm = _jsPlumb.getRenderMode();
        if (_jp.Endpoints[rm][t]) return new _jp.Endpoints[rm][t](p);
        if (!_jsPlumb.Defaults.DoNotThrowErrors)
          throw { msg: "jsPlumb: unknown endpoint type '" + t + "'" };
      };
      var endpointArgs = {
        _jsPlumb: this._jsPlumb.instance,
        cssClass: params$jscomp$0.cssClass,
        container: params$jscomp$0.container,
        tooltip: params$jscomp$0.tooltip,
        connectorTooltip: params$jscomp$0.connectorTooltip,
        endpoint: this,
      };
      if (_ju.isString(ep)) var endpoint = _e(ep, endpointArgs);
      else if (_ju.isArray(ep)) {
        endpointArgs = _ju.merge(ep[1], endpointArgs);
        endpoint = _e(ep[0], endpointArgs);
      } else endpoint = ep.clone();
      endpoint.clone = function () {
        if (_ju.isString(ep)) return _e(ep, endpointArgs);
        else if (_ju.isArray(ep)) {
          endpointArgs = _ju.merge(ep[1], endpointArgs);
          return _e(ep[0], endpointArgs);
        }
      }.bind(this);
      endpoint.typeId = typeId;
      return endpoint;
    };
    this.setEndpoint = function (ep, doNotRepaint) {
      var _ep = this.prepareEndpoint(ep);
      this.setPreparedEndpoint(_ep, true);
    };
    this.setPreparedEndpoint = function (ep, doNotRepaint) {
      if (this.endpoint != null) {
        this.endpoint.cleanup();
        this.endpoint.destroy();
      }
      this.endpoint = ep;
      this.type = this.endpoint.type;
      this.canvas = this.endpoint.canvas;
    };
    _jp.extend(this, params$jscomp$0, typeParameters);
    this.isSource = params$jscomp$0.isSource || false;
    this.isTemporarySource = params$jscomp$0.isTemporarySource || false;
    this.isTarget = params$jscomp$0.isTarget || false;
    this.connections = params$jscomp$0.connections || [];
    this.connectorPointerEvents = params$jscomp$0["connector-pointer-events"];
    this.scope = params$jscomp$0.scope || _jsPlumb.getDefaultScope();
    this.timestamp = null;
    this.reattachConnections =
      params$jscomp$0.reattach || _jsPlumb.Defaults.ReattachConnections;
    this.connectionsDetachable = _jsPlumb.Defaults.ConnectionsDetachable;
    if (
      params$jscomp$0.connectionsDetachable === false ||
      params$jscomp$0.detachable === false
    )
      this.connectionsDetachable = false;
    this.dragAllowedWhenFull = params$jscomp$0.dragAllowedWhenFull !== false;
    if (params$jscomp$0.onMaxConnections)
      this.bind("maxConnections", params$jscomp$0.onMaxConnections);
    this.addConnection = function (connection) {
      this.connections.push(connection);
      this[(this.connections.length > 0 ? "add" : "remove") + "Class"](
        _jsPlumb.endpointConnectedClass
      );
      this[(this.isFull() ? "add" : "remove") + "Class"](
        _jsPlumb.endpointFullClass
      );
    };
    this.detachFromConnection = function (connection, idx, doNotCleanup) {
      idx = idx == null ? this.connections.indexOf(connection) : idx;
      if (idx >= 0) {
        this.connections.splice(idx, 1);
        this[(this.connections.length > 0 ? "add" : "remove") + "Class"](
          _jsPlumb.endpointConnectedClass
        );
        this[(this.isFull() ? "add" : "remove") + "Class"](
          _jsPlumb.endpointFullClass
        );
      }
      if (!doNotCleanup && deleteOnEmpty && this.connections.length === 0)
        _jsPlumb.deleteObject({
          endpoint: this,
          fireEvent: false,
          deleteAttachedObjects: doNotCleanup !== true,
        });
    };
    this.deleteEveryConnection = function (params) {
      var c = this.connections.length;
      for (var i = 0; i < c; i++)
        _jsPlumb.deleteConnection(this.connections[0], params);
    };
    this.detachFrom = function (targetEndpoint, fireEvent, originalEvent) {
      var c = [];
      for (var i = 0; i < this.connections.length; i++)
        if (
          this.connections[i].endpoints[1] === targetEndpoint ||
          this.connections[i].endpoints[0] === targetEndpoint
        )
          c.push(this.connections[i]);
      var j = 0;
      for (var count = c.length; j < count; j++)
        _jsPlumb.deleteConnection(c[0]);
      return this;
    };
    this.getElement = function () {
      return this.element;
    };
    this.setElement = function (el) {
      var parentId = this._jsPlumb.instance.getId(el);
      var curId = this.elementId;
      _ju.removeWithFunction(
        params$jscomp$0.endpointsByElement[this.elementId],
        function (e) {
          return e.id === this.id;
        }.bind(this)
      );
      this.element = _jp.getElement(el);
      this.elementId = _jsPlumb.getId(this.element);
      _jsPlumb.router.rehomeEndpoint(this, curId, this.element);
      _jsPlumb.dragManager.endpointAdded(this.element);
      _ju.addToList(params$jscomp$0.endpointsByElement, parentId, this);
      return this;
    };
    this.makeInPlaceCopy = function () {
      var loc = this.anchor.getCurrentLocation({ element: this });
      var o = this.anchor.getOrientation(this);
      var acc = this.anchor.getCssClass();
      var inPlaceAnchor = {
        bind: function () {},
        compute: function () {
          return [loc[0], loc[1]];
        },
        getCurrentLocation: function () {
          return [loc[0], loc[1]];
        },
        getOrientation: function () {
          return o;
        },
        getCssClass: function () {
          return acc;
        },
      };
      return _newEndpoint({
        dropOptions: params$jscomp$0.dropOptions,
        anchor: inPlaceAnchor,
        source: this.element,
        paintStyle: this.getPaintStyle(),
        endpoint: params$jscomp$0.hideOnDrag ? "Blank" : this.endpoint,
        _transient: true,
        scope: this.scope,
        reference: this,
      });
    };
    this.connectorSelector = function () {
      return this.connections[0];
    };
    this.setStyle = this.setPaintStyle;
    this.paint = function (params) {
      params = params || {};
      var timestamp = params.timestamp;
      var recalc = !(params.recalc === false);
      if (!timestamp || this.timestamp !== timestamp) {
        var info = _jsPlumb.updateOffset({ elId: this.elementId, timestamp });
        var xy = params.offset ? params.offset.o : info.o;
        if (xy != null) {
          var ap = params.anchorPoint;
          var connectorPaintStyle = params.connectorPaintStyle;
          if (ap == null) {
            var wh = params.dimensions || info.s;
            var anchorParams = {
              xy: [xy.left, xy.top],
              wh,
              element: this,
              timestamp,
            };
            if (
              recalc &&
              this.anchor.isDynamic &&
              this.connections.length > 0
            ) {
              var c = findConnectionToUseForDynamicAnchor(
                this,
                params.elementWithPrecedence
              );
              var oIdx = c.endpoints[0] === this ? 1 : 0;
              var oId = oIdx === 0 ? c.sourceId : c.targetId;
              var oInfo = _jsPlumb.getCachedData(oId);
              var oOffset = oInfo.o;
              var oWH = oInfo.s;
              anchorParams.index = oIdx === 0 ? 1 : 0;
              anchorParams.connection = c;
              anchorParams.txy = [oOffset.left, oOffset.top];
              anchorParams.twh = oWH;
              anchorParams.tElement = c.endpoints[oIdx];
              anchorParams.tRotation = _jsPlumb.getRotation(oId);
            } else if (this.connections.length > 0)
              anchorParams.connection = this.connections[0];
            anchorParams.rotation = _jsPlumb.getRotation(this.elementId);
            ap = this.anchor.compute(anchorParams);
          }
          this.endpoint.compute(
            ap,
            this.anchor.getOrientation(this),
            this._jsPlumb.paintStyleInUse,
            connectorPaintStyle || this.paintStyleInUse
          );
          this.endpoint.paint(this._jsPlumb.paintStyleInUse, this.anchor);
          this.timestamp = timestamp;
          for (var i in this._jsPlumb.overlays)
            if (this._jsPlumb.overlays.hasOwnProperty(i)) {
              var o = this._jsPlumb.overlays[i];
              if (o.isVisible()) {
                this._jsPlumb.overlayPlacements[i] = o.draw(
                  this.endpoint,
                  this._jsPlumb.paintStyleInUse
                );
                o.paint(this._jsPlumb.overlayPlacements[i]);
              }
            }
        }
      }
    };
    this.getTypeDescriptor = function () {
      return "endpoint";
    };
    this.isVisible = function () {
      return this._jsPlumb.visible;
    };
    this.repaint = this.paint;
    var draggingInitialised = false;
    this.initDraggable = function () {
      if (!draggingInitialised && _jp.isDragSupported(this.element)) {
        var placeholderInfo = { id: null, element: null };
        var jpc = null;
        var existingJpc = false;
        var existingJpcParams = null;
        var _dragHandler = _makeConnectionDragHandler(
          this,
          placeholderInfo,
          _jsPlumb
        );
        var dragOptions = params$jscomp$0.dragOptions || {};
        var defaultOpts = {};
        var startEvent = _jp.dragEvents.start;
        var stopEvent = _jp.dragEvents.stop;
        var dragEvent = _jp.dragEvents.drag;
        var beforeStartEvent = _jp.dragEvents.beforeStart;
        var payload;
        var beforeStart = function (beforeStartParams) {
          payload = beforeStartParams.e.payload || {};
        };
        var start = function (startParams) {
          jpc = this.connectorSelector();
          var _continue = true;
          if (!this.isEnabled()) _continue = false;
          if (jpc == null && !this.isSource && !this.isTemporarySource)
            _continue = false;
          if (
            this.isSource &&
            this.isFull() &&
            !(jpc != null && this.dragAllowedWhenFull)
          )
            _continue = false;
          if (jpc != null && !jpc.isDetachable(this))
            if (this.isFull()) _continue = false;
            else jpc = null;
          var beforeDrag = _jsPlumb.checkCondition(
            jpc == null ? "beforeDrag" : "beforeStartDetach",
            {
              endpoint: this,
              source: this.element,
              sourceId: this.elementId,
              connection: jpc,
            }
          );
          if (beforeDrag === false) _continue = false;
          else if (typeof beforeDrag === "object")
            _jp.extend(beforeDrag, payload || {});
          else beforeDrag = payload || {};
          if (_continue === false) {
            if (_jsPlumb.stopDrag) _jsPlumb.stopDrag(this.canvas);
            _dragHandler.stopDrag();
            return false;
          }
          for (var i = 0; i < this.connections.length; i++)
            this.connections[i].setHover(false);
          this.addClass("endpointDrag");
          _jsPlumb.setConnectionBeingDragged(true);
          if (jpc && !this.isFull() && this.isSource) jpc = null;
          _jsPlumb.updateOffset({ elId: this.elementId });
          var ipco = this._jsPlumb.instance.getOffset(this.canvas);
          var canvasElement = this.canvas;
          var ips = this._jsPlumb.instance.getSize(this.canvas);
          _makeDraggablePlaceholder(placeholderInfo, _jsPlumb, ipco, ips);
          _jsPlumb.setAttributes(this.canvas, {
            dragId: placeholderInfo.id,
            elId: this.elementId,
          });
          var endpointToFloat = this.dragProxy || this.endpoint;
          if (this.dragProxy == null && this.connectionType != null) {
            var aae = this._jsPlumb.instance.deriveEndpointAndAnchorSpec(
              this.connectionType
            );
            if (aae.endpoints[1]) endpointToFloat = aae.endpoints[1];
          }
          var centerAnchor = this._jsPlumb.instance.makeAnchor("Center");
          centerAnchor.isFloating = true;
          this._jsPlumb.floatingEndpoint = _makeFloatingEndpoint(
            this.getPaintStyle(),
            centerAnchor,
            endpointToFloat,
            this.canvas,
            placeholderInfo.element,
            _jsPlumb,
            _newEndpoint,
            this.scope
          );
          var _savedAnchor = this._jsPlumb.floatingEndpoint.anchor;
          if (jpc == null) {
            this.setHover(false, false);
            jpc = _newConnection({
              sourceEndpoint: this,
              targetEndpoint: this._jsPlumb.floatingEndpoint,
              source: this.element,
              target: placeholderInfo.element,
              anchors: [this.anchor, this._jsPlumb.floatingEndpoint.anchor],
              paintStyle: params$jscomp$0.connectorStyle,
              hoverPaintStyle: params$jscomp$0.connectorHoverStyle,
              connector: params$jscomp$0.connector,
              overlays: params$jscomp$0.connectorOverlays,
              type: this.connectionType,
              cssClass: this.connectorClass,
              hoverClass: this.connectorHoverClass,
              scope: params$jscomp$0.scope,
              data: beforeDrag,
            });
            jpc.pending = true;
            jpc.addClass(_jsPlumb.draggingClass);
            this._jsPlumb.floatingEndpoint.addClass(_jsPlumb.draggingClass);
            this._jsPlumb.floatingEndpoint.anchor = _savedAnchor;
            _jsPlumb.fire("connectionDrag", jpc);
            _jsPlumb.router.newConnection(jpc);
          } else {
            existingJpc = true;
            jpc.setHover(false);
            var anchorIdx = jpc.endpoints[0].id === this.id ? 0 : 1;
            this.detachFromConnection(jpc, null, true);
            var dragScope = _jsPlumb.getDragScope(canvasElement);
            _jsPlumb.setAttribute(this.canvas, "originalScope", dragScope);
            _jsPlumb.fire("connectionDrag", jpc);
            if (anchorIdx === 0) {
              existingJpcParams = [
                jpc.source,
                jpc.sourceId,
                canvasElement,
                dragScope,
              ];
              _jsPlumb.router.sourceOrTargetChanged(
                jpc.endpoints[anchorIdx].elementId,
                placeholderInfo.id,
                jpc,
                placeholderInfo.element,
                0
              );
            } else {
              existingJpcParams = [
                jpc.target,
                jpc.targetId,
                canvasElement,
                dragScope,
              ];
              _jsPlumb.router.sourceOrTargetChanged(
                jpc.endpoints[anchorIdx].elementId,
                placeholderInfo.id,
                jpc,
                placeholderInfo.element,
                1
              );
            }
            jpc.suspendedEndpoint = jpc.endpoints[anchorIdx];
            jpc.suspendedElement = jpc.endpoints[anchorIdx].getElement();
            jpc.suspendedElementId = jpc.endpoints[anchorIdx].elementId;
            jpc.suspendedElementType = anchorIdx === 0 ? "source" : "target";
            jpc.suspendedEndpoint.setHover(false);
            this._jsPlumb.floatingEndpoint.referenceEndpoint =
              jpc.suspendedEndpoint;
            jpc.endpoints[anchorIdx] = this._jsPlumb.floatingEndpoint;
            jpc.addClass(_jsPlumb.draggingClass);
            this._jsPlumb.floatingEndpoint.addClass(_jsPlumb.draggingClass);
          }
          _jsPlumb.registerFloatingConnection(
            placeholderInfo,
            jpc,
            this._jsPlumb.floatingEndpoint
          );
          _jsPlumb.currentlyDragging = true;
        }.bind(this);
        var stop = function () {
          _jsPlumb.setConnectionBeingDragged(false);
          if (jpc && jpc.endpoints != null) {
            var originalEvent = _jsPlumb.getDropEvent(arguments);
            var idx = _jsPlumb.getFloatingAnchorIndex(jpc);
            jpc.endpoints[idx === 0 ? 1 : 0].anchor.locked = false;
            jpc.removeClass(_jsPlumb.draggingClass);
            if (
              this._jsPlumb &&
              (jpc.deleteConnectionNow ||
                jpc.endpoints[idx] === this._jsPlumb.floatingEndpoint)
            )
              if (existingJpc && jpc.suspendedEndpoint) {
                if (idx === 0) {
                  jpc.floatingElement = jpc.source;
                  jpc.floatingId = jpc.sourceId;
                  jpc.floatingEndpoint = jpc.endpoints[0];
                  jpc.floatingIndex = 0;
                  jpc.source = existingJpcParams[0];
                  jpc.sourceId = existingJpcParams[1];
                } else {
                  jpc.floatingElement = jpc.target;
                  jpc.floatingId = jpc.targetId;
                  jpc.floatingEndpoint = jpc.endpoints[1];
                  jpc.floatingIndex = 1;
                  jpc.target = existingJpcParams[0];
                  jpc.targetId = existingJpcParams[1];
                }
                var fe = this._jsPlumb.floatingEndpoint;
                _jsPlumb.setDragScope(
                  existingJpcParams[2],
                  existingJpcParams[3]
                );
                jpc.endpoints[idx] = jpc.suspendedEndpoint;
                if (
                  jpc.isReattach() ||
                  jpc._forceReattach ||
                  jpc._forceDetach ||
                  !_jsPlumb.deleteConnection(jpc, { originalEvent })
                ) {
                  jpc.setHover(false);
                  jpc._forceDetach = null;
                  jpc._forceReattach = null;
                  this._jsPlumb.floatingEndpoint.detachFromConnection(jpc);
                  jpc.suspendedEndpoint.addConnection(jpc);
                  if (idx === 1)
                    _jsPlumb.router.sourceOrTargetChanged(
                      jpc.floatingId,
                      jpc.targetId,
                      jpc,
                      jpc.target,
                      idx
                    );
                  else
                    _jsPlumb.router.sourceOrTargetChanged(
                      jpc.floatingId,
                      jpc.sourceId,
                      jpc,
                      jpc.source,
                      idx
                    );
                  _jsPlumb.repaint(existingJpcParams[1]);
                } else _jsPlumb.deleteObject({ endpoint: fe });
              }
            if (this.deleteAfterDragStop)
              _jsPlumb.deleteObject({ endpoint: this });
            else if (this._jsPlumb) this.paint({ recalc: false });
            _jsPlumb.fire("connectionDragStop", jpc, originalEvent);
            if (jpc.pending)
              _jsPlumb.fire("connectionAborted", jpc, originalEvent);
            _jsPlumb.currentlyDragging = false;
            jpc.suspendedElement = null;
            jpc.suspendedEndpoint = null;
            jpc = null;
          }
          if (placeholderInfo && placeholderInfo.element)
            _jsPlumb.remove(placeholderInfo.element, false, false);
          if (inPlaceCopy) _jsPlumb.deleteObject({ endpoint: inPlaceCopy });
          if (this._jsPlumb) {
            this.canvas.style.visibility = "visible";
            this.anchor.locked = false;
            this._jsPlumb.floatingEndpoint = null;
          }
        }.bind(this);
        dragOptions = _jp.extend(defaultOpts, dragOptions);
        dragOptions.scope = this.scope || dragOptions.scope;
        dragOptions[beforeStartEvent] = _ju.wrap(
          dragOptions[beforeStartEvent],
          beforeStart,
          false
        );
        dragOptions[startEvent] = _ju.wrap(
          dragOptions[startEvent],
          start,
          false
        );
        dragOptions[dragEvent] = _ju.wrap(
          dragOptions[dragEvent],
          _dragHandler.drag
        );
        dragOptions[stopEvent] = _ju.wrap(dragOptions[stopEvent], stop);
        dragOptions.multipleDrop = false;
        dragOptions.canDrag = function () {
          return (
            this.isSource ||
            this.isTemporarySource ||
            (this.connections.length > 0 &&
              this.connectionsDetachable !== false)
          );
        }.bind(this);
        _jsPlumb.initDraggable(this.canvas, dragOptions, "internal");
        this.canvas._jsPlumbRelatedElement = this.element;
        draggingInitialised = true;
      }
    };
    var ep$jscomp$0 =
      params$jscomp$0.endpoint ||
      this._jsPlumb.instance.Defaults.Endpoint ||
      _jp.Defaults.Endpoint;
    this.setEndpoint(ep$jscomp$0, true);
    var anchorParamsToUse = params$jscomp$0.anchor
      ? params$jscomp$0.anchor
      : params$jscomp$0.anchors
      ? params$jscomp$0.anchors
      : _jsPlumb.Defaults.Anchor || "Top";
    this.setAnchor(anchorParamsToUse, true);
    var type = ["default", params$jscomp$0.type || ""].join(" ");
    this.addType(type, params$jscomp$0.data, true);
    this.canvas = this.endpoint.canvas;
    this.canvas._jsPlumb = this;
    this.initDraggable();
    var _initDropTarget = function (
      canvas,
      isTransient,
      endpoint,
      referenceEndpoint
    ) {
      if (_jp.isDropSupported(this.element)) {
        var dropOptions =
          params$jscomp$0.dropOptions ||
          _jsPlumb.Defaults.DropOptions ||
          _jp.Defaults.DropOptions;
        dropOptions = _jp.extend({}, dropOptions);
        dropOptions.scope = dropOptions.scope || this.scope;
        var dropEvent = _jp.dragEvents.drop;
        var overEvent = _jp.dragEvents.over;
        var outEvent = _jp.dragEvents.out;
        var _ep = this;
        var drop = _jsPlumb.EndpointDropHandler({
          getEndpoint: function () {
            return _ep;
          },
          jsPlumb: _jsPlumb,
          enabled: function () {
            return endpoint != null ? endpoint.isEnabled() : true;
          },
          isFull: function () {
            return endpoint.isFull();
          },
          element: this.element,
          elementId: this.elementId,
          isSource: this.isSource,
          isTarget: this.isTarget,
          addClass: function (clazz) {
            _ep.addClass(clazz);
          },
          removeClass: function (clazz) {
            _ep.removeClass(clazz);
          },
          isDropAllowed: function () {
            return _ep.isDropAllowed.apply(_ep, arguments);
          },
          reference: referenceEndpoint,
          isRedrop: function (jpc, dhParams) {
            return (
              jpc.suspendedEndpoint &&
              dhParams.reference &&
              jpc.suspendedEndpoint.id === dhParams.reference.id
            );
          },
        });
        dropOptions[dropEvent] = _ju.wrap(dropOptions[dropEvent], drop, true);
        dropOptions[overEvent] = _ju.wrap(
          dropOptions[overEvent],
          function () {
            var draggable = _jp.getDragObject(arguments);
            var id = _jsPlumb.getAttribute(_jp.getElement(draggable), "dragId");
            var _jpc = _jsPlumb.getFloatingConnectionFor(id);
            if (_jpc != null) {
              var idx = _jsPlumb.getFloatingAnchorIndex(_jpc);
              var _cont =
                (this.isTarget && idx !== 0) ||
                (_jpc.suspendedEndpoint &&
                  this.referenceEndpoint &&
                  this.referenceEndpoint.id === _jpc.suspendedEndpoint.id);
              if (_cont) {
                var bb = _jsPlumb.checkCondition("checkDropAllowed", {
                  sourceEndpoint: _jpc.endpoints[idx],
                  targetEndpoint: this,
                  connection: _jpc,
                });
                this[(bb ? "add" : "remove") + "Class"](
                  _jsPlumb.endpointDropAllowedClass
                );
                this[(bb ? "remove" : "add") + "Class"](
                  _jsPlumb.endpointDropForbiddenClass
                );
                _jpc.endpoints[idx].anchor.over(this.anchor, this);
              }
            }
          }.bind(this)
        );
        dropOptions[outEvent] = _ju.wrap(
          dropOptions[outEvent],
          function () {
            var draggable = _jp.getDragObject(arguments);
            var id =
              draggable == null
                ? null
                : _jsPlumb.getAttribute(_jp.getElement(draggable), "dragId");
            var _jpc = id ? _jsPlumb.getFloatingConnectionFor(id) : null;
            if (_jpc != null) {
              var idx = _jsPlumb.getFloatingAnchorIndex(_jpc);
              var _cont =
                (this.isTarget && idx !== 0) ||
                (_jpc.suspendedEndpoint &&
                  this.referenceEndpoint &&
                  this.referenceEndpoint.id === _jpc.suspendedEndpoint.id);
              if (_cont) {
                this.removeClass(_jsPlumb.endpointDropAllowedClass);
                this.removeClass(_jsPlumb.endpointDropForbiddenClass);
                _jpc.endpoints[idx].anchor.out();
              }
            }
          }.bind(this)
        );
        _jsPlumb.initDroppable(canvas, dropOptions, "internal", isTransient);
      }
    }.bind(this);
    if (!this.anchor.isFloating)
      _initDropTarget(
        this.canvas,
        !(params$jscomp$0._transient || this.anchor.isFloating),
        this,
        params$jscomp$0.reference
      );
    return this;
  };
  _ju.extend(_jp.Endpoint, _jp.OverlayCapableJsPlumbUIComponent, {
    setVisible: function (v, doNotChangeConnections, doNotNotifyOtherEndpoint) {
      this._jsPlumb.visible = v;
      if (this.canvas) this.canvas.style.display = v ? "block" : "none";
      this[v ? "showOverlays" : "hideOverlays"]();
      if (!doNotChangeConnections)
        for (var i = 0; i < this.connections.length; i++) {
          this.connections[i].setVisible(v);
          if (!doNotNotifyOtherEndpoint) {
            var oIdx = this === this.connections[i].endpoints[0] ? 1 : 0;
            if (this.connections[i].endpoints[oIdx].connections.length === 1)
              this.connections[i].endpoints[oIdx].setVisible(v, true, true);
          }
        }
    },
    getAttachedElements: function () {
      return this.connections;
    },
    applyType: function (t, doNotRepaint) {
      this.setPaintStyle(t.endpointStyle || t.paintStyle, doNotRepaint);
      this.setHoverPaintStyle(
        t.endpointHoverStyle || t.hoverPaintStyle,
        doNotRepaint
      );
      if (t.maxConnections != null)
        this._jsPlumb.maxConnections = t.maxConnections;
      if (t.scope) this.scope = t.scope;
      _jp.extend(this, t, typeParameters);
      if (t.cssClass != null && this.canvas)
        this._jsPlumb.instance.addClass(this.canvas, t.cssClass);
      _jp.OverlayCapableJsPlumbUIComponent.applyType(this, t);
    },
    isEnabled: function () {
      return this._jsPlumb.enabled;
    },
    setEnabled: function (e) {
      this._jsPlumb.enabled = e;
    },
    cleanup: function () {
      var anchorClass =
        this._jsPlumb.instance.endpointAnchorClassPrefix +
        (this._jsPlumb.currentAnchorClass
          ? "-" + this._jsPlumb.currentAnchorClass
          : "");
      _jp.removeClass(this.element, anchorClass);
      this.anchor = null;
      this.endpoint.cleanup(true);
      this.endpoint.destroy();
      this.endpoint = null;
      this._jsPlumb.instance.destroyDraggable(this.canvas, "internal");
      this._jsPlumb.instance.destroyDroppable(this.canvas, "internal");
    },
    setHover: function (h) {
      if (
        this.endpoint &&
        this._jsPlumb &&
        !this._jsPlumb.instance.isConnectionBeingDragged()
      )
        this.endpoint.setHover(h);
    },
    isFull: function () {
      return this._jsPlumb.maxConnections === 0
        ? true
        : !(
            this.isFloating() ||
            this._jsPlumb.maxConnections < 0 ||
            this.connections.length < this._jsPlumb.maxConnections
          );
    },
    isFloating: function () {
      return this.anchor != null && this.anchor.isFloating;
    },
    isConnectedTo: function (endpoint) {
      var found = false;
      if (endpoint)
        for (var i = 0; i < this.connections.length; i++)
          if (
            this.connections[i].endpoints[1] === endpoint ||
            this.connections[i].endpoints[0] === endpoint
          ) {
            found = true;
            break;
          }
      return found;
    },
    getConnectionCost: function () {
      return this._jsPlumb.connectionCost;
    },
    setConnectionCost: function (c) {
      this._jsPlumb.connectionCost = c;
    },
    areConnectionsDirected: function () {
      return this._jsPlumb.connectionsDirected;
    },
    setConnectionsDirected: function (b) {
      this._jsPlumb.connectionsDirected = b;
    },
    setElementId: function (_elId) {
      this.elementId = _elId;
      this.anchor.elementId = _elId;
    },
    setReferenceElement: function (_el) {
      this.element = _jp.getElement(_el);
    },
    setDragAllowedWhenFull: function (allowed) {
      this.dragAllowedWhenFull = allowed;
    },
    equals: function (endpoint) {
      return this.anchor.equals(endpoint.anchor);
    },
    getUuid: function () {
      return this._jsPlumb.uuid;
    },
    computeAnchor: function (params) {
      return this.anchor.compute(params);
    },
  });
  root.jsPlumbInstance.prototype.EndpointDropHandler = function (dhParams) {
    return function (e) {
      var _jsPlumb = dhParams.jsPlumb;
      dhParams.removeClass(_jsPlumb.endpointDropAllowedClass);
      dhParams.removeClass(_jsPlumb.endpointDropForbiddenClass);
      var originalEvent = _jsPlumb.getDropEvent(arguments);
      var draggable = _jsPlumb.getDragObject(arguments);
      var id = _jsPlumb.getAttribute(draggable, "dragId");
      var elId = _jsPlumb.getAttribute(draggable, "elId");
      var scope = _jsPlumb.getAttribute(draggable, "originalScope");
      var jpc = _jsPlumb.getFloatingConnectionFor(id);
      if (jpc == null) return;
      var existingConnection = jpc.suspendedEndpoint != null;
      if (existingConnection && jpc.suspendedEndpoint._jsPlumb == null) return;
      var _ep = dhParams.getEndpoint(jpc);
      if (_ep == null) return;
      if (dhParams.isRedrop(jpc, dhParams)) {
        jpc._forceReattach = true;
        jpc.setHover(false);
        if (dhParams.maybeCleanup) dhParams.maybeCleanup(_ep);
        return;
      }
      var idx = _jsPlumb.getFloatingAnchorIndex(jpc);
      if (
        (idx === 0 && !dhParams.isSource) ||
        (idx === 1 && !dhParams.isTarget)
      ) {
        if (dhParams.maybeCleanup) dhParams.maybeCleanup(_ep);
        return;
      }
      if (dhParams.onDrop) dhParams.onDrop(jpc);
      if (scope) _jsPlumb.setDragScope(draggable, scope);
      var isFull = dhParams.isFull(e);
      if (isFull)
        _ep.fire(
          "maxConnections",
          {
            endpoint: this,
            connection: jpc,
            maxConnections: _ep._jsPlumb.maxConnections,
          },
          originalEvent
        );
      if (!isFull && dhParams.enabled()) {
        var _doContinue = true;
        if (idx === 0) {
          jpc.floatingElement = jpc.source;
          jpc.floatingId = jpc.sourceId;
          jpc.floatingEndpoint = jpc.endpoints[0];
          jpc.floatingIndex = 0;
          jpc.source = dhParams.element;
          jpc.sourceId = _jsPlumb.getId(dhParams.element);
        } else {
          jpc.floatingElement = jpc.target;
          jpc.floatingId = jpc.targetId;
          jpc.floatingEndpoint = jpc.endpoints[1];
          jpc.floatingIndex = 1;
          jpc.target = dhParams.element;
          jpc.targetId = _jsPlumb.getId(dhParams.element);
        }
        if (existingConnection && jpc.suspendedEndpoint.id !== _ep.id)
          if (
            !jpc.isDetachAllowed(jpc) ||
            !jpc.endpoints[idx].isDetachAllowed(jpc) ||
            !jpc.suspendedEndpoint.isDetachAllowed(jpc) ||
            !_jsPlumb.checkCondition("beforeDetach", jpc)
          )
            _doContinue = false;
        var continueFunction = function (optionalData) {
          jpc.endpoints[idx].detachFromConnection(jpc);
          if (jpc.suspendedEndpoint)
            jpc.suspendedEndpoint.detachFromConnection(jpc);
          jpc.endpoints[idx] = _ep;
          _ep.addConnection(jpc);
          var params = _ep.getParameters();
          for (var aParam in params) jpc.setParameter(aParam, params[aParam]);
          if (!existingConnection) {
            if (params.draggable)
              _jsPlumb.initDraggable(
                this.element,
                dhParams.dragOptions,
                "internal",
                _jsPlumb
              );
          } else {
            var suspendedElementId = jpc.suspendedEndpoint.elementId;
            _jsPlumb.fireMoveEvent(
              {
                index: idx,
                originalSourceId: idx === 0 ? suspendedElementId : jpc.sourceId,
                newSourceId: idx === 0 ? _ep.elementId : jpc.sourceId,
                originalTargetId: idx === 1 ? suspendedElementId : jpc.targetId,
                newTargetId: idx === 1 ? _ep.elementId : jpc.targetId,
                originalSourceEndpoint:
                  idx === 0 ? jpc.suspendedEndpoint : jpc.endpoints[0],
                newSourceEndpoint: idx === 0 ? _ep : jpc.endpoints[0],
                originalTargetEndpoint:
                  idx === 1 ? jpc.suspendedEndpoint : jpc.endpoints[1],
                newTargetEndpoint: idx === 1 ? _ep : jpc.endpoints[1],
                connection: jpc,
              },
              originalEvent
            );
          }
          if (idx === 1)
            _jsPlumb.router.sourceOrTargetChanged(
              jpc.floatingId,
              jpc.targetId,
              jpc,
              jpc.target,
              1
            );
          else
            _jsPlumb.router.sourceOrTargetChanged(
              jpc.floatingId,
              jpc.sourceId,
              jpc,
              jpc.source,
              0
            );
          if (jpc.endpoints[0].finalEndpoint) {
            var _toDelete = jpc.endpoints[0];
            _toDelete.detachFromConnection(jpc);
            jpc.endpoints[0] = jpc.endpoints[0].finalEndpoint;
            jpc.endpoints[0].addConnection(jpc);
          }
          if (_ju.isObject(optionalData)) jpc.mergeData(optionalData);
          _jsPlumb.finaliseConnection(jpc, null, originalEvent, false);
          jpc.setHover(false);
          _jsPlumb.revalidate(jpc.endpoints[0].element);
        }.bind(this);
        var dontContinueFunction = function () {
          if (jpc.suspendedEndpoint) {
            jpc.endpoints[idx] = jpc.suspendedEndpoint;
            jpc.setHover(false);
            jpc._forceDetach = true;
            if (idx === 0) {
              jpc.source = jpc.suspendedEndpoint.element;
              jpc.sourceId = jpc.suspendedEndpoint.elementId;
            } else {
              jpc.target = jpc.suspendedEndpoint.element;
              jpc.targetId = jpc.suspendedEndpoint.elementId;
            }
            jpc.suspendedEndpoint.addConnection(jpc);
            if (idx === 1)
              _jsPlumb.router.sourceOrTargetChanged(
                jpc.floatingId,
                jpc.targetId,
                jpc,
                jpc.target,
                1
              );
            else
              _jsPlumb.router.sourceOrTargetChanged(
                jpc.floatingId,
                jpc.sourceId,
                jpc,
                jpc.source,
                0
              );
            _jsPlumb.repaint(jpc.sourceId);
            jpc._forceDetach = false;
          }
        };
        _doContinue =
          _doContinue &&
          dhParams.isDropAllowed(
            jpc.sourceId,
            jpc.targetId,
            jpc.scope,
            jpc,
            _ep
          );
        if (_doContinue) {
          continueFunction(_doContinue);
          return true;
        } else dontContinueFunction();
      }
      if (dhParams.maybeCleanup) dhParams.maybeCleanup(_ep);
      _jsPlumb.currentlyDragging = false;
    };
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var makeConnector = function (
    _jsPlumb,
    renderMode,
    connectorName,
    connectorArgs,
    forComponent
  ) {
    _jp.Connectors[renderMode] = _jp.Connectors[renderMode] || {};
    if (_jp.Connectors[renderMode][connectorName] == null) {
      if (_jp.Connectors[connectorName] == null)
        if (!_jsPlumb.Defaults.DoNotThrowErrors)
          throw new TypeError(
            "jsPlumb: unknown connector type '" + connectorName + "'"
          );
        else return null;
      _jp.Connectors[renderMode][connectorName] = function () {
        _jp.Connectors[connectorName].apply(this, arguments);
        _jp.ConnectorRenderers[renderMode].apply(this, arguments);
      };
      _ju.extend(_jp.Connectors[renderMode][connectorName], [
        _jp.Connectors[connectorName],
        _jp.ConnectorRenderers[renderMode],
      ]);
    }
    return new _jp.Connectors[renderMode][connectorName](
      connectorArgs,
      forComponent
    );
  };
  var _makeAnchor = function (anchorParams, elementId, _jsPlumb) {
    return anchorParams
      ? _jsPlumb.makeAnchor(anchorParams, elementId, _jsPlumb)
      : null;
  };
  var _updateConnectedClass = function (conn, element, _jsPlumb, remove) {
    if (element != null) {
      element._jsPlumbConnections = element._jsPlumbConnections || {};
      if (remove) delete element._jsPlumbConnections[conn.id];
      else element._jsPlumbConnections[conn.id] = true;
      if (_ju.isEmpty(element._jsPlumbConnections))
        _jsPlumb.removeClass(element, _jsPlumb.connectedClass);
      else _jsPlumb.addClass(element, _jsPlumb.connectedClass);
    }
  };
  _jp.Connection = function (params) {
    var _newEndpoint = params.newEndpoint;
    this.id = params.id;
    this.connector = null;
    this.idPrefix = "_jsplumb_c_";
    this.defaultLabelLocation = 0.5;
    this.defaultOverlayKeys = ["Overlays", "ConnectionOverlays"];
    this.previousConnection = params.previousConnection;
    this.source = _jp.getElement(params.source);
    this.target = _jp.getElement(params.target);
    _jp.OverlayCapableJsPlumbUIComponent.apply(this, arguments);
    if (params.sourceEndpoint) {
      this.source = params.sourceEndpoint.getElement();
      this.sourceId = params.sourceEndpoint.elementId;
    } else this.sourceId = this._jsPlumb.instance.getId(this.source);
    if (params.targetEndpoint) {
      this.target = params.targetEndpoint.getElement();
      this.targetId = params.targetEndpoint.elementId;
    } else this.targetId = this._jsPlumb.instance.getId(this.target);
    this.scope = params.scope;
    this.endpoints = [];
    this.endpointStyles = [];
    var _jsPlumb = this._jsPlumb.instance;
    _jsPlumb.manage(this.sourceId, this.source);
    _jsPlumb.manage(this.targetId, this.target);
    this._jsPlumb.visible = true;
    this._jsPlumb.params = {
      cssClass: params.cssClass,
      container: params.container,
      "pointer-events": params["pointer-events"],
      editorParams: params.editorParams,
      overlays: params.overlays,
    };
    this._jsPlumb.lastPaintedAt = null;
    this.bind(
      "mouseover",
      function () {
        this.setHover(true);
      }.bind(this)
    );
    this.bind(
      "mouseout",
      function () {
        this.setHover(false);
      }.bind(this)
    );
    this.makeEndpoint = function (isSource, el, elId, ep, definition) {
      elId = elId || this._jsPlumb.instance.getId(el);
      return this.prepareEndpoint(
        _jsPlumb,
        _newEndpoint,
        this,
        ep,
        isSource ? 0 : 1,
        params,
        el,
        elId,
        definition
      );
    };
    if (params.type)
      params.endpoints =
        params.endpoints ||
        this._jsPlumb.instance.deriveEndpointAndAnchorSpec(params.type)
          .endpoints;
    var eS = this.makeEndpoint(
      true,
      this.source,
      this.sourceId,
      params.sourceEndpoint
    );
    var eT = this.makeEndpoint(
      false,
      this.target,
      this.targetId,
      params.targetEndpoint
    );
    if (eS) _ju.addToList(params.endpointsByElement, this.sourceId, eS);
    if (eT) _ju.addToList(params.endpointsByElement, this.targetId, eT);
    if (!this.scope) this.scope = this.endpoints[0].scope;
    if (params.deleteEndpointsOnEmpty != null) {
      this.endpoints[0].setDeleteOnEmpty(params.deleteEndpointsOnEmpty);
      this.endpoints[1].setDeleteOnEmpty(params.deleteEndpointsOnEmpty);
    }
    var _detachable = _jsPlumb.Defaults.ConnectionsDetachable;
    if (params.detachable === false) _detachable = false;
    if (this.endpoints[0].connectionsDetachable === false) _detachable = false;
    if (this.endpoints[1].connectionsDetachable === false) _detachable = false;
    var _reattach =
      params.reattach ||
      this.endpoints[0].reattachConnections ||
      this.endpoints[1].reattachConnections ||
      _jsPlumb.Defaults.ReattachConnections;
    this.appendToDefaultType({
      detachable: _detachable,
      reattach: _reattach,
      paintStyle:
        this.endpoints[0].connectorStyle ||
        this.endpoints[1].connectorStyle ||
        params.paintStyle ||
        _jsPlumb.Defaults.PaintStyle ||
        _jp.Defaults.PaintStyle,
      hoverPaintStyle:
        this.endpoints[0].connectorHoverStyle ||
        this.endpoints[1].connectorHoverStyle ||
        params.hoverPaintStyle ||
        _jsPlumb.Defaults.HoverPaintStyle ||
        _jp.Defaults.HoverPaintStyle,
    });
    var _suspendedAt = _jsPlumb.getSuspendedAt();
    if (!_jsPlumb.isSuspendDrawing()) {
      var myInfo = _jsPlumb.getCachedData(this.sourceId);
      var myOffset = myInfo.o;
      var myWH = myInfo.s;
      var otherInfo = _jsPlumb.getCachedData(this.targetId);
      var otherOffset = otherInfo.o;
      var otherWH = otherInfo.s;
      var initialTimestamp = _suspendedAt || jsPlumbUtil.uuid();
      var anchorLoc = this.endpoints[0].anchor.compute({
        xy: [myOffset.left, myOffset.top],
        wh: myWH,
        element: this.endpoints[0],
        elementId: this.endpoints[0].elementId,
        txy: [otherOffset.left, otherOffset.top],
        twh: otherWH,
        tElement: this.endpoints[1],
        timestamp: initialTimestamp,
        rotation: _jsPlumb.getRotation(this.endpoints[0].elementId),
      });
      this.endpoints[0].paint({ anchorLoc, timestamp: initialTimestamp });
      anchorLoc = this.endpoints[1].anchor.compute({
        xy: [otherOffset.left, otherOffset.top],
        wh: otherWH,
        element: this.endpoints[1],
        elementId: this.endpoints[1].elementId,
        txy: [myOffset.left, myOffset.top],
        twh: myWH,
        tElement: this.endpoints[0],
        timestamp: initialTimestamp,
        rotation: _jsPlumb.getRotation(this.endpoints[1].elementId),
      });
      this.endpoints[1].paint({ anchorLoc, timestamp: initialTimestamp });
    }
    this.getTypeDescriptor = function () {
      return "connection";
    };
    this.getAttachedElements = function () {
      return this.endpoints;
    };
    this.isDetachable = function (ep) {
      return this._jsPlumb.detachable === false
        ? false
        : ep != null
        ? ep.connectionsDetachable === true
        : this._jsPlumb.detachable === true;
    };
    this.setDetachable = function (detachable) {
      this._jsPlumb.detachable = detachable === true;
    };
    this.isReattach = function () {
      return (
        this._jsPlumb.reattach === true ||
        this.endpoints[0].reattachConnections === true ||
        this.endpoints[1].reattachConnections === true
      );
    };
    this.setReattach = function (reattach) {
      this._jsPlumb.reattach = reattach === true;
    };
    this._jsPlumb.cost = params.cost || this.endpoints[0].getConnectionCost();
    this._jsPlumb.directed = params.directed;
    if (params.directed == null)
      this._jsPlumb.directed = this.endpoints[0].areConnectionsDirected();
    var _p = _jp.extend({}, this.endpoints[1].getParameters());
    _jp.extend(_p, this.endpoints[0].getParameters());
    _jp.extend(_p, this.getParameters());
    this.setParameters(_p);
    this.setConnector(
      this.endpoints[0].connector ||
        this.endpoints[1].connector ||
        params.connector ||
        _jsPlumb.Defaults.Connector ||
        _jp.Defaults.Connector,
      true
    );
    var data =
      params.data == null || !_ju.isObject(params.data) ? {} : params.data;
    this.getData = function () {
      return data;
    };
    this.setData = function (d) {
      data = d || {};
    };
    this.mergeData = function (d) {
      data = _jp.extend(data, d);
    };
    var _types = [
      "default",
      this.endpoints[0].connectionType,
      this.endpoints[1].connectionType,
      params.type,
    ].join(" ");
    if (/[^\s]/.test(_types)) this.addType(_types, params.data, true);
    this.updateConnectedClass();
  };
  _ju.extend(_jp.Connection, _jp.OverlayCapableJsPlumbUIComponent, {
    applyType: function (t, doNotRepaint, typeMap) {
      var _connector = null;
      if (t.connector != null) {
        _connector = this.getCachedTypeItem("connector", typeMap.connector);
        if (_connector == null) {
          _connector = this.prepareConnector(t.connector, typeMap.connector);
          this.cacheTypeItem("connector", _connector, typeMap.connector);
        }
        this.setPreparedConnector(_connector);
      }
      if (t.detachable != null) this.setDetachable(t.detachable);
      if (t.reattach != null) this.setReattach(t.reattach);
      if (t.scope) this.scope = t.scope;
      if (t.cssClass != null && this.canvas)
        this._jsPlumb.instance.addClass(this.canvas, t.cssClass);
      var _anchors = null;
      if (t.anchor) {
        _anchors = this.getCachedTypeItem("anchors", typeMap.anchor);
        if (_anchors == null) {
          _anchors = [
            this._jsPlumb.instance.makeAnchor(t.anchor),
            this._jsPlumb.instance.makeAnchor(t.anchor),
          ];
          this.cacheTypeItem("anchors", _anchors, typeMap.anchor);
        }
      } else if (t.anchors) {
        _anchors = this.getCachedTypeItem("anchors", typeMap.anchors);
        if (_anchors == null) {
          _anchors = [
            this._jsPlumb.instance.makeAnchor(t.anchors[0]),
            this._jsPlumb.instance.makeAnchor(t.anchors[1]),
          ];
          this.cacheTypeItem("anchors", _anchors, typeMap.anchors);
        }
      }
      if (_anchors != null) {
        this.endpoints[0].anchor = _anchors[0];
        this.endpoints[1].anchor = _anchors[1];
        if (this.endpoints[1].anchor.isDynamic)
          this._jsPlumb.instance.repaint(this.endpoints[1].elementId);
      }
      _jp.OverlayCapableJsPlumbUIComponent.applyType(this, t);
    },
    addClass: function (c, informEndpoints) {
      if (informEndpoints) {
        this.endpoints[0].addClass(c);
        this.endpoints[1].addClass(c);
        if (this.suspendedEndpoint) this.suspendedEndpoint.addClass(c);
      }
      if (this.connector) this.connector.addClass(c);
    },
    removeClass: function (c, informEndpoints) {
      if (informEndpoints) {
        this.endpoints[0].removeClass(c);
        this.endpoints[1].removeClass(c);
        if (this.suspendedEndpoint) this.suspendedEndpoint.removeClass(c);
      }
      if (this.connector) this.connector.removeClass(c);
    },
    isVisible: function () {
      return this._jsPlumb.visible;
    },
    setVisible: function (v) {
      this._jsPlumb.visible = v;
      if (this.connector) this.connector.setVisible(v);
      this.repaint();
    },
    cleanup: function () {
      this.updateConnectedClass(true);
      this.endpoints = null;
      this.source = null;
      this.target = null;
      if (this.connector != null) {
        this.connector.cleanup(true);
        this.connector.destroy(true);
      }
      this.connector = null;
    },
    updateConnectedClass: function (remove) {
      if (this._jsPlumb) {
        _updateConnectedClass(
          this,
          this.source,
          this._jsPlumb.instance,
          remove
        );
        _updateConnectedClass(
          this,
          this.target,
          this._jsPlumb.instance,
          remove
        );
      }
    },
    setHover: function (state) {
      if (
        this.connector &&
        this._jsPlumb &&
        !this._jsPlumb.instance.isConnectionBeingDragged()
      ) {
        this.connector.setHover(state);
        root.jsPlumb[state ? "addClass" : "removeClass"](
          this.source,
          this._jsPlumb.instance.hoverSourceClass
        );
        root.jsPlumb[state ? "addClass" : "removeClass"](
          this.target,
          this._jsPlumb.instance.hoverTargetClass
        );
      }
    },
    getUuids: function () {
      return [this.endpoints[0].getUuid(), this.endpoints[1].getUuid()];
    },
    getCost: function () {
      return this._jsPlumb ? this._jsPlumb.cost : -Infinity;
    },
    setCost: function (c) {
      this._jsPlumb.cost = c;
    },
    isDirected: function () {
      return this._jsPlumb.directed;
    },
    getConnector: function () {
      return this.connector;
    },
    prepareConnector: function (connectorSpec, typeId) {
      var connectorArgs = {
        _jsPlumb: this._jsPlumb.instance,
        cssClass: this._jsPlumb.params.cssClass,
        container: this._jsPlumb.params.container,
        "pointer-events": this._jsPlumb.params["pointer-events"],
      };
      var renderMode = this._jsPlumb.instance.getRenderMode();
      if (_ju.isString(connectorSpec))
        var connector = makeConnector(
          this._jsPlumb.instance,
          renderMode,
          connectorSpec,
          connectorArgs,
          this
        );
      else if (_ju.isArray(connectorSpec))
        if (connectorSpec.length === 1)
          connector = makeConnector(
            this._jsPlumb.instance,
            renderMode,
            connectorSpec[0],
            connectorArgs,
            this
          );
        else
          connector = makeConnector(
            this._jsPlumb.instance,
            renderMode,
            connectorSpec[0],
            _ju.merge(connectorSpec[1], connectorArgs),
            this
          );
      if (typeId != null) connector.typeId = typeId;
      return connector;
    },
    setPreparedConnector: function (
      connector,
      doNotRepaint,
      doNotChangeListenerComponent,
      typeId
    ) {
      if (this.connector !== connector) {
        var previousClasses = "";
        if (this.connector != null) {
          var previous = this.connector;
          previousClasses = previous.getClass();
          this.connector.cleanup();
          this.connector.destroy();
        }
        this.connector = connector;
        if (typeId) this.cacheTypeItem("connector", connector, typeId);
        this.canvas = this.connector.canvas;
        this.bgCanvas = this.connector.bgCanvas;
        this.connector.reattach(this._jsPlumb.instance);
        this.addClass(previousClasses);
        if (this.canvas) this.canvas._jsPlumb = this;
        if (this.bgCanvas) this.bgCanvas._jsPlumb = this;
        if (previous != null) {
          var o = this.getOverlays();
          for (var i = 0; i < o.length; i++)
            if (o[i].transfer) o[i].transfer(this.connector);
        }
        if (!doNotChangeListenerComponent)
          this.setListenerComponent(this.connector);
        if (!doNotRepaint) this.repaint();
      }
    },
    setConnector: function (
      connectorSpec,
      doNotRepaint,
      doNotChangeListenerComponent,
      typeId
    ) {
      var connector = this.prepareConnector(connectorSpec, typeId);
      this.setPreparedConnector(
        connector,
        doNotRepaint,
        doNotChangeListenerComponent,
        typeId
      );
    },
    paint: function (params) {
      if (!this._jsPlumb.instance.isSuspendDrawing() && this._jsPlumb.visible) {
        params = params || {};
        var timestamp = params.timestamp;
        var swap = false;
        var tId = swap ? this.sourceId : this.targetId;
        var sId = swap ? this.targetId : this.sourceId;
        var tIdx = swap ? 0 : 1;
        var sIdx = swap ? 1 : 0;
        if (timestamp == null || timestamp !== this._jsPlumb.lastPaintedAt) {
          var sourceInfo = this._jsPlumb.instance.updateOffset({ elId: sId }).o;
          var targetInfo = this._jsPlumb.instance.updateOffset({ elId: tId }).o;
          var sE = this.endpoints[sIdx];
          var tE = this.endpoints[tIdx];
          var sAnchorP = sE.anchor.getCurrentLocation({
            xy: [sourceInfo.left, sourceInfo.top],
            wh: [sourceInfo.width, sourceInfo.height],
            element: sE,
            timestamp,
            rotation: this._jsPlumb.instance.getRotation(this.sourceId),
          });
          var tAnchorP = tE.anchor.getCurrentLocation({
            xy: [targetInfo.left, targetInfo.top],
            wh: [targetInfo.width, targetInfo.height],
            element: tE,
            timestamp,
            rotation: this._jsPlumb.instance.getRotation(this.targetId),
          });
          this.connector.resetBounds();
          this.connector.compute({
            sourcePos: sAnchorP,
            targetPos: tAnchorP,
            sourceOrientation: sE.anchor.getOrientation(sE),
            targetOrientation: tE.anchor.getOrientation(tE),
            sourceEndpoint: this.endpoints[sIdx],
            targetEndpoint: this.endpoints[tIdx],
            "stroke-width": this._jsPlumb.paintStyleInUse.strokeWidth,
            sourceInfo,
            targetInfo,
          });
          var overlayExtents = {
            minX: Infinity,
            minY: Infinity,
            maxX: -Infinity,
            maxY: -Infinity,
          };
          for (var i in this._jsPlumb.overlays)
            if (this._jsPlumb.overlays.hasOwnProperty(i)) {
              var o = this._jsPlumb.overlays[i];
              if (o.isVisible()) {
                this._jsPlumb.overlayPlacements[i] = o.draw(
                  this.connector,
                  this._jsPlumb.paintStyleInUse,
                  this.getAbsoluteOverlayPosition(o)
                );
                overlayExtents.minX = Math.min(
                  overlayExtents.minX,
                  this._jsPlumb.overlayPlacements[i].minX
                );
                overlayExtents.maxX = Math.max(
                  overlayExtents.maxX,
                  this._jsPlumb.overlayPlacements[i].maxX
                );
                overlayExtents.minY = Math.min(
                  overlayExtents.minY,
                  this._jsPlumb.overlayPlacements[i].minY
                );
                overlayExtents.maxY = Math.max(
                  overlayExtents.maxY,
                  this._jsPlumb.overlayPlacements[i].maxY
                );
              }
            }
          var lineWidth =
            parseFloat(this._jsPlumb.paintStyleInUse.strokeWidth || 1) / 2;
          var outlineWidth = parseFloat(
            this._jsPlumb.paintStyleInUse.strokeWidth || 0
          );
          var extents = {
            xmin: Math.min(
              this.connector.bounds.minX - (lineWidth + outlineWidth),
              overlayExtents.minX
            ),
            ymin: Math.min(
              this.connector.bounds.minY - (lineWidth + outlineWidth),
              overlayExtents.minY
            ),
            xmax: Math.max(
              this.connector.bounds.maxX + (lineWidth + outlineWidth),
              overlayExtents.maxX
            ),
            ymax: Math.max(
              this.connector.bounds.maxY + (lineWidth + outlineWidth),
              overlayExtents.maxY
            ),
          };
          this.connector.paintExtents = extents;
          this.connector.paint(this._jsPlumb.paintStyleInUse, null, extents);
          for (var j in this._jsPlumb.overlays)
            if (this._jsPlumb.overlays.hasOwnProperty(j)) {
              var p = this._jsPlumb.overlays[j];
              if (p.isVisible())
                p.paint(this._jsPlumb.overlayPlacements[j], extents);
            }
        }
        this._jsPlumb.lastPaintedAt = timestamp;
      }
    },
    repaint: function (params) {
      var p = jsPlumb.extend(params || {}, {});
      p.elId = this.sourceId;
      this.paint(p);
    },
    prepareEndpoint: function (
      _jsPlumb,
      _newEndpoint,
      conn,
      existing,
      index,
      params,
      element,
      elementId,
      definition
    ) {
      if (existing) {
        conn.endpoints[index] = existing;
        existing.addConnection(conn);
      } else {
        if (!params.endpoints) params.endpoints = [null, null];
        var ep =
          definition ||
          params.endpoints[index] ||
          params.endpoint ||
          _jsPlumb.Defaults.Endpoints[index] ||
          _jp.Defaults.Endpoints[index] ||
          _jsPlumb.Defaults.Endpoint ||
          _jp.Defaults.Endpoint;
        if (!params.endpointStyles) params.endpointStyles = [null, null];
        if (!params.endpointHoverStyles)
          params.endpointHoverStyles = [null, null];
        var es =
          params.endpointStyles[index] ||
          params.endpointStyle ||
          _jsPlumb.Defaults.EndpointStyles[index] ||
          _jp.Defaults.EndpointStyles[index] ||
          _jsPlumb.Defaults.EndpointStyle ||
          _jp.Defaults.EndpointStyle;
        if (es.fill == null && params.paintStyle != null)
          es.fill = params.paintStyle.stroke;
        if (es.outlineStroke == null && params.paintStyle != null)
          es.outlineStroke = params.paintStyle.outlineStroke;
        if (es.outlineWidth == null && params.paintStyle != null)
          es.outlineWidth = params.paintStyle.outlineWidth;
        var ehs =
          params.endpointHoverStyles[index] ||
          params.endpointHoverStyle ||
          _jsPlumb.Defaults.EndpointHoverStyles[index] ||
          _jp.Defaults.EndpointHoverStyles[index] ||
          _jsPlumb.Defaults.EndpointHoverStyle ||
          _jp.Defaults.EndpointHoverStyle;
        if (params.hoverPaintStyle != null) {
          if (ehs == null) ehs = {};
          if (ehs.fill == null) ehs.fill = params.hoverPaintStyle.stroke;
        }
        var a = params.anchors
          ? params.anchors[index]
          : params.anchor
          ? params.anchor
          : _makeAnchor(
              _jsPlumb.Defaults.Anchors[index],
              elementId,
              _jsPlumb
            ) ||
            _makeAnchor(_jp.Defaults.Anchors[index], elementId, _jsPlumb) ||
            _makeAnchor(_jsPlumb.Defaults.Anchor, elementId, _jsPlumb) ||
            _makeAnchor(_jp.Defaults.Anchor, elementId, _jsPlumb);
        var u = params.uuids ? params.uuids[index] : null;
        var e = _newEndpoint({
          paintStyle: es,
          hoverPaintStyle: ehs,
          endpoint: ep,
          connections: [conn],
          uuid: u,
          anchor: a,
          source: element,
          scope: params.scope,
          reattach: params.reattach || _jsPlumb.Defaults.ReattachConnections,
          detachable:
            params.detachable || _jsPlumb.Defaults.ConnectionsDetachable,
        });
        if (existing == null) e.setDeleteOnEmpty(true);
        conn.endpoints[index] = e;
        if (params.drawEndpoints === false) e.setVisible(false, true, true);
      }
      return e;
    },
    replaceEndpoint: function (idx, endpointDef) {
      var current = this.endpoints[idx];
      var elId = current.elementId;
      var ebe = this._jsPlumb.instance.getEndpoints(elId);
      var _idx = ebe.indexOf(current);
      var _new = this.makeEndpoint(
        idx === 0,
        current.element,
        elId,
        null,
        endpointDef
      );
      this.endpoints[idx] = _new;
      ebe.splice(_idx, 1, _new);
      this._jsPlumb.instance.deleteObject({
        endpoint: current,
        deleteAttachedObjects: false,
      });
      this._jsPlumb.instance.fire("endpointReplaced", {
        previous: current,
        current: _new,
      });
      this._jsPlumb.instance.router.sourceOrTargetChanged(
        this.endpoints[1].elementId,
        this.endpoints[1].elementId,
        this,
        this.endpoints[1].element,
        1
      );
    },
  });
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _ju = root.jsPlumbUtil;
  var _jp = root.jsPlumb;
  _jp.AnchorManager = function (params$jscomp$0) {
    var _amEndpoints = {};
    var continuousAnchorLocations = {};
    var continuousAnchorOrientations = {};
    var connectionsByElementId = {};
    var self = this;
    var anchorLists = {};
    var jsPlumbInstance = params$jscomp$0.jsPlumbInstance;
    var floatingConnections = {};
    var placeAnchorsOnLine = function (
      desc,
      elementDimensions,
      elementPosition,
      connections,
      horizontal,
      otherMultiplier,
      reverse,
      rotation
    ) {
      var a = [];
      var step =
        elementDimensions[horizontal ? 0 : 1] / (connections.length + 1);
      for (var i = 0; i < connections.length; i++) {
        var val = (i + 1) * step;
        var other = otherMultiplier * elementDimensions[horizontal ? 1 : 0];
        if (reverse) val = elementDimensions[horizontal ? 0 : 1] - val;
        var dx = horizontal ? val : other;
        var x = elementPosition.left + dx;
        var xp = dx / elementDimensions[0];
        var dy = horizontal ? other : val;
        var y = elementPosition.top + dy;
        var yp = dy / elementDimensions[1];
        if (rotation !== 0) {
          var rotated = jsPlumbUtil.rotatePoint(
            [x, y],
            [elementPosition.centerx, elementPosition.centery],
            rotation
          );
          x = rotated[0];
          y = rotated[1];
        }
        a.push([x, y, xp, yp, connections[i][1], connections[i][2]]);
      }
      return a;
    };
    var rightAndBottomSort = function (a, b) {
      return b[0][0] - a[0][0];
    };
    var leftAndTopSort = function (a, b) {
      var p1 = a[0][0] < 0 ? -Math.PI - a[0][0] : Math.PI - a[0][0];
      var p2 = b[0][0] < 0 ? -Math.PI - b[0][0] : Math.PI - b[0][0];
      return p1 - p2;
    };
    var edgeSortFunctions = {
      top: leftAndTopSort,
      right: rightAndBottomSort,
      bottom: rightAndBottomSort,
      left: leftAndTopSort,
    };
    var _sortHelper = function (_array, _fn) {
      return _array.sort(_fn);
    };
    var placeAnchors = function (elementId, _anchorLists) {
      var cd = jsPlumbInstance.getCachedData(elementId);
      var sS = cd.s;
      var sO = cd.o;
      var placeSomeAnchors = function (
        desc,
        elementDimensions,
        elementPosition,
        unsortedConnections,
        isHorizontal,
        otherMultiplier,
        orientation
      ) {
        if (unsortedConnections.length > 0) {
          var sc = _sortHelper(unsortedConnections, edgeSortFunctions[desc]);
          var reverse = desc === "right" || desc === "top";
          var rotation = jsPlumbInstance.getRotation(elementId);
          var anchors = placeAnchorsOnLine(
            desc,
            elementDimensions,
            elementPosition,
            sc,
            isHorizontal,
            otherMultiplier,
            reverse,
            rotation
          );
          var _setAnchorLocation = function (endpoint, anchorPos) {
            continuousAnchorLocations[endpoint.id] = [
              anchorPos[0],
              anchorPos[1],
              anchorPos[2],
              anchorPos[3],
            ];
            continuousAnchorOrientations[endpoint.id] = orientation;
          };
          for (var i = 0; i < anchors.length; i++) {
            var c = anchors[i][4];
            var weAreSource = c.endpoints[0].elementId === elementId;
            var weAreTarget = c.endpoints[1].elementId === elementId;
            if (weAreSource) _setAnchorLocation(c.endpoints[0], anchors[i]);
            if (weAreTarget) _setAnchorLocation(c.endpoints[1], anchors[i]);
          }
        }
      };
      placeSomeAnchors("bottom", sS, sO, _anchorLists.bottom, true, 1, [0, 1]);
      placeSomeAnchors("top", sS, sO, _anchorLists.top, true, 0, [0, -1]);
      placeSomeAnchors("left", sS, sO, _anchorLists.left, false, 0, [-1, 0]);
      placeSomeAnchors("right", sS, sO, _anchorLists.right, false, 1, [1, 0]);
    };
    this.reset = function () {
      _amEndpoints = {};
      connectionsByElementId = {};
      anchorLists = {};
    };
    this.addFloatingConnection = function (key, conn) {
      floatingConnections[key] = conn;
    };
    this.newConnection = function (conn) {
      var sourceId = conn.sourceId;
      var targetId = conn.targetId;
      var ep = conn.endpoints;
      var doRegisterTarget = true;
      var registerConnection = function (
        otherIndex,
        otherEndpoint,
        otherAnchor,
        elId,
        c
      ) {
        if (sourceId === targetId && otherAnchor.isContinuous) {
          conn._jsPlumb.instance.removeElement(ep[1].canvas);
          doRegisterTarget = false;
        }
        _ju.addToList(connectionsByElementId, elId, [
          c,
          otherEndpoint,
          otherAnchor.constructor === _jp.DynamicAnchor,
        ]);
      };
      registerConnection(0, ep[0], ep[0].anchor, targetId, conn);
      if (doRegisterTarget)
        registerConnection(1, ep[1], ep[1].anchor, sourceId, conn);
    };
    var removeEndpointFromAnchorLists = function (endpoint) {
      (function (list, eId) {
        if (list) {
          var f = function (e) {
            return e[4] === eId;
          };
          _ju.removeWithFunction(list.top, f);
          _ju.removeWithFunction(list.left, f);
          _ju.removeWithFunction(list.bottom, f);
          _ju.removeWithFunction(list.right, f);
        }
      })(anchorLists[endpoint.elementId], endpoint.id);
    };
    this.connectionDetached = function (connInfo, doNotRedraw) {
      var connection = connInfo.connection || connInfo;
      var sourceId = connInfo.sourceId;
      var targetId = connInfo.targetId;
      var ep = connection.endpoints;
      var removeConnection = function (
        otherIndex,
        otherEndpoint,
        otherAnchor,
        elId,
        c
      ) {
        _ju.removeWithFunction(connectionsByElementId[elId], function (_c) {
          return _c[0].id === c.id;
        });
      };
      removeConnection(1, ep[1], ep[1].anchor, sourceId, connection);
      removeConnection(0, ep[0], ep[0].anchor, targetId, connection);
      if (connection.floatingId) {
        removeConnection(
          connection.floatingIndex,
          connection.floatingEndpoint,
          connection.floatingEndpoint.anchor,
          connection.floatingId,
          connection
        );
        removeEndpointFromAnchorLists(connection.floatingEndpoint);
      }
      removeEndpointFromAnchorLists(connection.endpoints[0]);
      removeEndpointFromAnchorLists(connection.endpoints[1]);
      if (!doNotRedraw) {
        self.redraw(connection.sourceId);
        if (connection.targetId !== connection.sourceId)
          self.redraw(connection.targetId);
      }
    };
    this.addEndpoint = function (endpoint, elementId) {
      _ju.addToList(_amEndpoints, elementId, endpoint);
    };
    this.changeId = function (oldId, newId) {
      connectionsByElementId[newId] = connectionsByElementId[oldId];
      _amEndpoints[newId] = _amEndpoints[oldId];
      delete connectionsByElementId[oldId];
      delete _amEndpoints[oldId];
    };
    this.getConnectionsFor = function (elementId) {
      return connectionsByElementId[elementId] || [];
    };
    this.getEndpointsFor = function (elementId) {
      return _amEndpoints[elementId] || [];
    };
    this.deleteEndpoint = function (endpoint) {
      _ju.removeWithFunction(_amEndpoints[endpoint.elementId], function (e) {
        return e.id === endpoint.id;
      });
      removeEndpointFromAnchorLists(endpoint);
    };
    this.elementRemoved = function (elementId) {
      delete floatingConnections[elementId];
      delete _amEndpoints[elementId];
      _amEndpoints[elementId] = [];
    };
    var _updateAnchorList = function (
      lists,
      theta,
      order,
      conn,
      aBoolean,
      otherElId,
      idx,
      reverse,
      edgeId,
      elId,
      connsToPaint,
      endpointsToPaint
    ) {
      var exactIdx = -1;
      var firstMatchingElIdx = -1;
      var endpoint = conn.endpoints[idx];
      var endpointId = endpoint.id;
      var oIdx = [1, 0][idx];
      var values = [[theta, order], conn, aBoolean, otherElId, endpointId];
      var listToAddTo = lists[edgeId];
      var listToRemoveFrom = endpoint._continuousAnchorEdge
        ? lists[endpoint._continuousAnchorEdge]
        : null;
      var i;
      if (listToRemoveFrom) {
        var rIdx = _ju.findWithFunction(listToRemoveFrom, function (e) {
          return e[4] === endpointId;
        });
        if (rIdx !== -1) {
          listToRemoveFrom.splice(rIdx, 1);
          for (i = 0; i < listToRemoveFrom.length; i++) {
            var candidate = listToRemoveFrom[i][1];
            _ju.addWithFunction(connsToPaint, candidate, function (c) {
              return c.id === candidate.id;
            });
            _ju.addWithFunction(
              endpointsToPaint,
              listToRemoveFrom[i][1].endpoints[idx],
              function (e) {
                return e.id === candidate.endpoints[idx].id;
              }
            );
            _ju.addWithFunction(
              endpointsToPaint,
              listToRemoveFrom[i][1].endpoints[oIdx],
              function (e) {
                return e.id === candidate.endpoints[oIdx].id;
              }
            );
          }
        }
      }
      for (i = 0; i < listToAddTo.length; i++) {
        candidate = listToAddTo[i][1];
        if (
          params$jscomp$0.idx === 1 &&
          listToAddTo[i][3] === otherElId &&
          firstMatchingElIdx === -1
        )
          firstMatchingElIdx = i;
        _ju.addWithFunction(connsToPaint, candidate, function (c) {
          return c.id === candidate.id;
        });
        _ju.addWithFunction(
          endpointsToPaint,
          listToAddTo[i][1].endpoints[idx],
          function (e) {
            return e.id === candidate.endpoints[idx].id;
          }
        );
        _ju.addWithFunction(
          endpointsToPaint,
          listToAddTo[i][1].endpoints[oIdx],
          function (e) {
            return e.id === candidate.endpoints[oIdx].id;
          }
        );
      }
      if (exactIdx !== -1) listToAddTo[exactIdx] = values;
      else {
        var insertIdx = reverse
          ? firstMatchingElIdx !== -1
            ? firstMatchingElIdx
            : 0
          : listToAddTo.length;
        listToAddTo.splice(insertIdx, 0, values);
      }
      endpoint._continuousAnchorEdge = edgeId;
    };
    this.sourceOrTargetChanged = function (
      originalId,
      newId,
      connection,
      newElement,
      anchorIndex
    ) {
      if (anchorIndex === 0) {
        if (originalId !== newId) {
          connection.sourceId = newId;
          connection.source = newElement;
          _ju.removeWithFunction(
            connectionsByElementId[originalId],
            function (info) {
              return info[0].id === connection.id;
            }
          );
          var tIdx = _ju.findWithFunction(
            connectionsByElementId[connection.targetId],
            function (i) {
              return i[0].id === connection.id;
            }
          );
          if (tIdx > -1) {
            connectionsByElementId[connection.targetId][tIdx][0] = connection;
            connectionsByElementId[connection.targetId][tIdx][1] =
              connection.endpoints[0];
            connectionsByElementId[connection.targetId][tIdx][2] =
              connection.endpoints[0].anchor.constructor === _jp.DynamicAnchor;
          }
          _ju.addToList(connectionsByElementId, newId, [
            connection,
            connection.endpoints[1],
            connection.endpoints[1].anchor.constructor === _jp.DynamicAnchor,
          ]);
          if (connection.endpoints[1].anchor.isContinuous)
            if (connection.source === connection.target)
              connection._jsPlumb.instance.removeElement(
                connection.endpoints[1].canvas
              );
            else if (connection.endpoints[1].canvas.parentNode == null)
              connection._jsPlumb.instance.appendElement(
                connection.endpoints[1].canvas
              );
          connection.updateConnectedClass();
        }
      } else if (anchorIndex === 1) {
        var sourceElId = connection.endpoints[0].elementId;
        connection.target = newElement;
        connection.targetId = newId;
        var sIndex = _ju.findWithFunction(
          connectionsByElementId[sourceElId],
          function (i) {
            return i[0].id === connection.id;
          }
        );
        var tIndex = _ju.findWithFunction(
          connectionsByElementId[originalId],
          function (i) {
            return i[0].id === connection.id;
          }
        );
        if (sIndex !== -1) {
          connectionsByElementId[sourceElId][sIndex][0] = connection;
          connectionsByElementId[sourceElId][sIndex][1] =
            connection.endpoints[1];
          connectionsByElementId[sourceElId][sIndex][2] =
            connection.endpoints[1].anchor.constructor === _jp.DynamicAnchor;
        }
        if (tIndex > -1) {
          connectionsByElementId[originalId].splice(tIndex, 1);
          _ju.addToList(connectionsByElementId, newId, [
            connection,
            connection.endpoints[0],
            connection.endpoints[0].anchor.constructor === _jp.DynamicAnchor,
          ]);
        }
        connection.updateConnectedClass();
      }
    };
    this.rehomeEndpoint = function (ep, currentId, element) {
      var eps = _amEndpoints[currentId] || [];
      var elementId = jsPlumbInstance.getId(element);
      if (elementId !== currentId) {
        var idx = eps.indexOf(ep);
        if (idx > -1) {
          var _ep = eps.splice(idx, 1)[0];
          self.add(_ep, elementId);
        }
      }
      for (var i = 0; i < ep.connections.length; i++)
        if (ep.connections[i].sourceId === currentId)
          self.sourceOrTargetChanged(
            currentId,
            ep.elementId,
            ep.connections[i],
            ep.element,
            0
          );
        else if (ep.connections[i].targetId === currentId)
          self.sourceOrTargetChanged(
            currentId,
            ep.elementId,
            ep.connections[i],
            ep.element,
            1
          );
    };
    this.redraw = function (
      elementId,
      ui,
      timestamp,
      offsetToUI,
      clearEdits,
      doNotRecalcEndpoint
    ) {
      var connectionsToPaint = [];
      var endpointsToPaint = [];
      var anchorsToUpdate = [];
      if (!jsPlumbInstance.isSuspendDrawing()) {
        var ep = _amEndpoints[elementId] || [];
        var endpointConnections = connectionsByElementId[elementId] || [];
        timestamp = timestamp || jsPlumbUtil.uuid();
        offsetToUI = offsetToUI || { left: 0, top: 0 };
        if (ui)
          ui = {
            left: ui.left + offsetToUI.left,
            top: ui.top + offsetToUI.top,
          };
        var myOffset = jsPlumbInstance.updateOffset({
          elId: elementId,
          offset: ui,
          recalc: false,
          timestamp,
        });
        var orientationCache = {};
        for (var i = 0; i < endpointConnections.length; i++) {
          var conn = endpointConnections[i][0];
          var sourceId = conn.sourceId;
          var targetId = conn.targetId;
          var sourceContinuous = conn.endpoints[0].anchor.isContinuous;
          var targetContinuous = conn.endpoints[1].anchor.isContinuous;
          if (sourceContinuous || targetContinuous) {
            var oKey = sourceId + "_" + targetId;
            var o = orientationCache[oKey];
            var oIdx = conn.sourceId === elementId ? 1 : 0;
            var targetRotation = jsPlumbInstance.getRotation(targetId);
            var sourceRotation = jsPlumbInstance.getRotation(sourceId);
            if (sourceContinuous && !anchorLists[sourceId])
              anchorLists[sourceId] = {
                top: [],
                right: [],
                bottom: [],
                left: [],
              };
            if (targetContinuous && !anchorLists[targetId])
              anchorLists[targetId] = {
                top: [],
                right: [],
                bottom: [],
                left: [],
              };
            if (elementId !== targetId)
              jsPlumbInstance.updateOffset({ elId: targetId, timestamp });
            if (elementId !== sourceId)
              jsPlumbInstance.updateOffset({ elId: sourceId, timestamp });
            var td = jsPlumbInstance.getCachedData(targetId);
            var sd = jsPlumbInstance.getCachedData(sourceId);
            if (
              targetId === sourceId &&
              (sourceContinuous || targetContinuous)
            ) {
              _updateAnchorList(
                anchorLists[sourceId],
                -Math.PI / 2,
                0,
                conn,
                false,
                targetId,
                0,
                false,
                "top",
                sourceId,
                connectionsToPaint,
                endpointsToPaint
              );
              _updateAnchorList(
                anchorLists[targetId],
                -Math.PI / 2,
                0,
                conn,
                false,
                sourceId,
                1,
                false,
                "top",
                targetId,
                connectionsToPaint,
                endpointsToPaint
              );
            } else {
              if (!o) {
                o = this.calculateOrientation(
                  sourceId,
                  targetId,
                  sd.o,
                  td.o,
                  conn.endpoints[0].anchor,
                  conn.endpoints[1].anchor,
                  conn,
                  sourceRotation,
                  targetRotation
                );
                orientationCache[oKey] = o;
              }
              if (sourceContinuous)
                _updateAnchorList(
                  anchorLists[sourceId],
                  o.theta,
                  0,
                  conn,
                  false,
                  targetId,
                  0,
                  false,
                  o.a[0],
                  sourceId,
                  connectionsToPaint,
                  endpointsToPaint
                );
              if (targetContinuous)
                _updateAnchorList(
                  anchorLists[targetId],
                  o.theta2,
                  -1,
                  conn,
                  true,
                  sourceId,
                  1,
                  true,
                  o.a[1],
                  targetId,
                  connectionsToPaint,
                  endpointsToPaint
                );
            }
            if (sourceContinuous)
              _ju.addWithFunction(anchorsToUpdate, sourceId, function (a) {
                return a === sourceId;
              });
            if (targetContinuous)
              _ju.addWithFunction(anchorsToUpdate, targetId, function (a) {
                return a === targetId;
              });
            _ju.addWithFunction(connectionsToPaint, conn, function (c) {
              return c.id === conn.id;
            });
            if (
              (sourceContinuous && oIdx === 0) ||
              (targetContinuous && oIdx === 1)
            )
              _ju.addWithFunction(
                endpointsToPaint,
                conn.endpoints[oIdx],
                function (e) {
                  return e.id === conn.endpoints[oIdx].id;
                }
              );
          }
        }
        for (i = 0; i < ep.length; i++)
          if (ep[i].connections.length === 0 && ep[i].anchor.isContinuous) {
            if (!anchorLists[elementId])
              anchorLists[elementId] = {
                top: [],
                right: [],
                bottom: [],
                left: [],
              };
            _updateAnchorList(
              anchorLists[elementId],
              -Math.PI / 2,
              0,
              { endpoints: [ep[i], ep[i]], paint: function () {} },
              false,
              elementId,
              0,
              false,
              ep[i].anchor.getDefaultFace(),
              elementId,
              connectionsToPaint,
              endpointsToPaint
            );
            _ju.addWithFunction(anchorsToUpdate, elementId, function (a) {
              return a === elementId;
            });
          }
        for (i = 0; i < anchorsToUpdate.length; i++)
          placeAnchors(anchorsToUpdate[i], anchorLists[anchorsToUpdate[i]]);
        for (i = 0; i < ep.length; i++)
          ep[i].paint({
            timestamp,
            offset: myOffset,
            dimensions: myOffset.s,
            recalc: doNotRecalcEndpoint !== true,
          });
        for (i = 0; i < endpointsToPaint.length; i++) {
          var cd = jsPlumbInstance.getCachedData(endpointsToPaint[i].elementId);
          endpointsToPaint[i].paint({
            timestamp: null,
            offset: cd,
            dimensions: cd.s,
          });
        }
        for (i = 0; i < endpointConnections.length; i++) {
          var otherEndpoint = endpointConnections[i][1];
          if (otherEndpoint.anchor.constructor === _jp.DynamicAnchor) {
            otherEndpoint.paint({
              elementWithPrecedence: elementId,
              timestamp,
            });
            _ju.addWithFunction(
              connectionsToPaint,
              endpointConnections[i][0],
              function (c) {
                return c.id === endpointConnections[i][0].id;
              }
            );
            for (var k = 0; k < otherEndpoint.connections.length; k++)
              if (otherEndpoint.connections[k] !== endpointConnections[i][0])
                _ju.addWithFunction(
                  connectionsToPaint,
                  otherEndpoint.connections[k],
                  function (c) {
                    return c.id === otherEndpoint.connections[k].id;
                  }
                );
          } else
            _ju.addWithFunction(
              connectionsToPaint,
              endpointConnections[i][0],
              function (c) {
                return c.id === endpointConnections[i][0].id;
              }
            );
        }
        var fc = floatingConnections[elementId];
        if (fc) fc.paint({ timestamp, recalc: false, elId: elementId });
        for (i = 0; i < connectionsToPaint.length; i++)
          connectionsToPaint[i].paint({
            elId: elementId,
            timestamp: null,
            recalc: false,
            clearEdits,
          });
      }
      return { c: connectionsToPaint, e: endpointsToPaint };
    };
    var ContinuousAnchor = function (anchorParams) {
      _ju.EventGenerator.apply(this);
      this.type = "Continuous";
      this.isDynamic = true;
      this.isContinuous = true;
      var faces = anchorParams.faces || ["top", "right", "bottom", "left"];
      var clockwise = !(anchorParams.clockwise === false);
      var availableFaces = {};
      var opposites = {
        top: "bottom",
        right: "left",
        left: "right",
        bottom: "top",
      };
      var clockwiseOptions = {
        top: "right",
        right: "bottom",
        left: "top",
        bottom: "left",
      };
      var antiClockwiseOptions = {
        top: "left",
        right: "top",
        left: "bottom",
        bottom: "right",
      };
      var secondBest = clockwise ? clockwiseOptions : antiClockwiseOptions;
      var lastChoice = clockwise ? antiClockwiseOptions : clockwiseOptions;
      var cssClass = anchorParams.cssClass || "";
      var _currentFace = null;
      var _lockedFace = null;
      var X_AXIS_FACES = ["left", "right"];
      var Y_AXIS_FACES = ["top", "bottom"];
      var _lockedAxis = null;
      for (var i = 0; i < faces.length; i++) availableFaces[faces[i]] = true;
      this.getDefaultFace = function () {
        return faces.length === 0 ? "top" : faces[0];
      };
      this.isRelocatable = function () {
        return true;
      };
      this.isSnapOnRelocate = function () {
        return true;
      };
      this.verifyEdge = function (edge) {
        if (availableFaces[edge]) return edge;
        else if (availableFaces[opposites[edge]]) return opposites[edge];
        else if (availableFaces[secondBest[edge]]) return secondBest[edge];
        else if (availableFaces[lastChoice[edge]]) return lastChoice[edge];
        return edge;
      };
      this.isEdgeSupported = function (edge) {
        return _lockedAxis == null
          ? _lockedFace == null
            ? availableFaces[edge] === true
            : _lockedFace === edge
          : _lockedAxis.indexOf(edge) !== -1;
      };
      this.setCurrentFace = function (face, overrideLock) {
        _currentFace = face;
        if (overrideLock && _lockedFace != null) _lockedFace = _currentFace;
      };
      this.getCurrentFace = function () {
        return _currentFace;
      };
      this.getSupportedFaces = function () {
        var af = [];
        for (var k in availableFaces) if (availableFaces[k]) af.push(k);
        return af;
      };
      this.lock = function () {
        _lockedFace = _currentFace;
      };
      this.unlock = function () {
        _lockedFace = null;
      };
      this.isLocked = function () {
        return _lockedFace != null;
      };
      this.lockCurrentAxis = function () {
        if (_currentFace != null)
          _lockedAxis =
            _currentFace === "left" || _currentFace === "right"
              ? X_AXIS_FACES
              : Y_AXIS_FACES;
      };
      this.unlockCurrentAxis = function () {
        _lockedAxis = null;
      };
      this.compute = function (params) {
        return continuousAnchorLocations[params.element.id] || [0, 0];
      };
      this.getCurrentLocation = function (params) {
        return continuousAnchorLocations[params.element.id] || [0, 0];
      };
      this.getOrientation = function (endpoint) {
        return continuousAnchorOrientations[endpoint.id] || [0, 0];
      };
      this.getCssClass = function () {
        return cssClass;
      };
    };
    jsPlumbInstance.continuousAnchorFactory = {
      get: function (params) {
        return new ContinuousAnchor(params);
      },
      clear: function (elementId) {
        delete continuousAnchorLocations[elementId];
      },
    };
  };
  _jp.AnchorManager.prototype.calculateOrientation = function (
    sourceId,
    targetId,
    sd,
    td,
    sourceAnchor,
    targetAnchor,
    connection,
    sourceRotation,
    targetRotation
  ) {
    var Orientation = {
      HORIZONTAL: "horizontal",
      VERTICAL: "vertical",
      DIAGONAL: "diagonal",
      IDENTITY: "identity",
    };
    var axes = ["left", "top", "right", "bottom"];
    if (sourceId === targetId)
      return { orientation: Orientation.IDENTITY, a: ["top", "top"] };
    var theta = Math.atan2(td.centery - sd.centery, td.centerx - sd.centerx);
    var theta2 = Math.atan2(sd.centery - td.centery, sd.centerx - td.centerx);
    var candidates = [];
    var midpoints = {};
    (function (types, dim) {
      for (var i = 0; i < types.length; i++) {
        midpoints[types[i]] = {
          left: [dim[i][0].left, dim[i][0].centery],
          right: [dim[i][0].right, dim[i][0].centery],
          top: [dim[i][0].centerx, dim[i][0].top],
          bottom: [dim[i][0].centerx, dim[i][0].bottom],
        };
        if (dim[i][1] !== 0)
          for (var axis in midpoints[types[i]])
            midpoints[types[i]][axis] = jsPlumbUtil.rotatePoint(
              midpoints[types[i]][axis],
              [dim[i][0].centerx, dim[i][0].centery],
              dim[i][1]
            );
      }
    })(
      ["source", "target"],
      [
        [sd, sourceRotation],
        [td, targetRotation],
      ]
    );
    for (var sf = 0; sf < axes.length; sf++)
      for (var tf = 0; tf < axes.length; tf++)
        candidates.push({
          source: axes[sf],
          target: axes[tf],
          dist: Biltong.lineLength(
            midpoints.source[axes[sf]],
            midpoints.target[axes[tf]]
          ),
        });
    candidates.sort(function (a, b) {
      return a.dist < b.dist ? -1 : a.dist > b.dist ? 1 : 0;
    });
    var sourceEdge = candidates[0].source;
    var targetEdge = candidates[0].target;
    for (var i$jscomp$0 = 0; i$jscomp$0 < candidates.length; i$jscomp$0++) {
      if (sourceAnchor.isContinuous && sourceAnchor.locked)
        sourceEdge = sourceAnchor.getCurrentFace();
      else if (
        !sourceAnchor.isContinuous ||
        sourceAnchor.isEdgeSupported(candidates[i$jscomp$0].source)
      )
        sourceEdge = candidates[i$jscomp$0].source;
      else sourceEdge = null;
      if (targetAnchor.isContinuous && targetAnchor.locked)
        targetEdge = targetAnchor.getCurrentFace();
      else if (
        !targetAnchor.isContinuous ||
        targetAnchor.isEdgeSupported(candidates[i$jscomp$0].target)
      )
        targetEdge = candidates[i$jscomp$0].target;
      else targetEdge = null;
      if (sourceEdge != null && targetEdge != null) break;
    }
    if (sourceAnchor.isContinuous) sourceAnchor.setCurrentFace(sourceEdge);
    if (targetAnchor.isContinuous) targetAnchor.setCurrentFace(targetEdge);
    return { a: [sourceEdge, targetEdge], theta, theta2 };
  };
  _jp.Anchor = function (params$jscomp$0) {
    this.x = params$jscomp$0.x || 0;
    this.y = params$jscomp$0.y || 0;
    this.elementId = params$jscomp$0.elementId;
    this.cssClass = params$jscomp$0.cssClass || "";
    this.orientation = params$jscomp$0.orientation || [0, 0];
    this.lastReturnValue = null;
    this.offsets = params$jscomp$0.offsets || [0, 0];
    this.timestamp = null;
    this._unrotatedOrientation = [this.orientation[0], this.orientation[1]];
    this.relocatable = params$jscomp$0.relocatable !== false;
    this.snapOnRelocate = params$jscomp$0.snapOnRelocate !== false;
    this.locked = false;
    _ju.EventGenerator.apply(this);
    this.compute = function (params) {
      var xy = params.xy;
      var wh = params.wh;
      var timestamp = params.timestamp;
      if (timestamp && timestamp === this.timestamp)
        return this.lastReturnValue;
      var candidate = [
        xy[0] + this.x * wh[0] + this.offsets[0],
        xy[1] + this.y * wh[1] + this.offsets[1],
        this.x,
        this.y,
      ];
      var rotation = params.rotation;
      if (rotation != null && rotation !== 0) {
        var c2 = jsPlumbUtil.rotatePoint(
          candidate,
          [xy[0] + wh[0] / 2, xy[1] + wh[1] / 2],
          rotation
        );
        this.orientation[0] = Math.round(
          this._unrotatedOrientation[0] * c2[2] -
            this._unrotatedOrientation[1] * c2[3]
        );
        this.orientation[1] = Math.round(
          this._unrotatedOrientation[1] * c2[2] +
            this._unrotatedOrientation[0] * c2[3]
        );
        this.lastReturnValue = [c2[0], c2[1], this.x, this.y];
      } else {
        this.orientation[0] = this._unrotatedOrientation[0];
        this.orientation[1] = this._unrotatedOrientation[1];
        this.lastReturnValue = candidate;
      }
      this.timestamp = timestamp;
      return this.lastReturnValue;
    };
    this.getCurrentLocation = function (params) {
      params = params || {};
      return this.lastReturnValue == null ||
        (params.timestamp != null && this.timestamp !== params.timestamp)
        ? this.compute(params)
        : this.lastReturnValue;
    };
    this.setPosition = function (x, y, ox, oy, overrideLock) {
      if (!this.locked || overrideLock) {
        this.x = x;
        this.y = y;
        this.orientation = [ox, oy];
        this.lastReturnValue = null;
      }
    };
  };
  _ju.extend(_jp.Anchor, _ju.EventGenerator, {
    equals: function (anchor) {
      if (!anchor) return false;
      var ao = anchor.getOrientation();
      var o = this.getOrientation();
      return (
        this.x === anchor.x &&
        this.y === anchor.y &&
        this.offsets[0] === anchor.offsets[0] &&
        this.offsets[1] === anchor.offsets[1] &&
        o[0] === ao[0] &&
        o[1] === ao[1]
      );
    },
    getOrientation: function () {
      return this.orientation;
    },
    getCssClass: function () {
      return this.cssClass;
    },
  });
  _jp.FloatingAnchor = function (params$jscomp$0) {
    _jp.Anchor.apply(this, arguments);
    var ref = params$jscomp$0.reference;
    var refCanvas = params$jscomp$0.referenceCanvas;
    var size = _jp.getSize(refCanvas);
    var xDir = 0;
    var yDir = 0;
    var orientation = null;
    var _lastResult = null;
    this.orientation = null;
    this.x = 0;
    this.y = 0;
    this.isFloating = true;
    this.compute = function (params) {
      var xy = params.xy;
      var result = [xy[0] + size[0] / 2, xy[1] + size[1] / 2];
      _lastResult = result;
      return result;
    };
    this.getOrientation = function (_endpoint) {
      if (orientation) return orientation;
      else {
        var o = ref.getOrientation(_endpoint);
        return [Math.abs(o[0]) * xDir * -1, Math.abs(o[1]) * yDir * -1];
      }
    };
    this.over = function (anchor, endpoint) {
      orientation = anchor.getOrientation(endpoint);
    };
    this.out = function () {
      orientation = null;
    };
    this.getCurrentLocation = function (params) {
      return _lastResult == null ? this.compute(params) : _lastResult;
    };
  };
  _ju.extend(_jp.FloatingAnchor, _jp.Anchor);
  var _convertAnchor = function (anchor, jsPlumbInstance, elementId) {
    return anchor.constructor === _jp.Anchor
      ? anchor
      : jsPlumbInstance.makeAnchor(anchor, elementId, jsPlumbInstance);
  };
  _jp.DynamicAnchor = function (params$jscomp$0) {
    _jp.Anchor.apply(this, arguments);
    this.isDynamic = true;
    this.anchors = [];
    this.elementId = params$jscomp$0.elementId;
    this.jsPlumbInstance = params$jscomp$0.jsPlumbInstance;
    for (
      var i$jscomp$0 = 0;
      i$jscomp$0 < params$jscomp$0.anchors.length;
      i$jscomp$0++
    )
      this.anchors[i$jscomp$0] = _convertAnchor(
        params$jscomp$0.anchors[i$jscomp$0],
        this.jsPlumbInstance,
        this.elementId
      );
    this.getAnchors = function () {
      return this.anchors;
    };
    var _curAnchor = this.anchors.length > 0 ? this.anchors[0] : null;
    var _lastAnchor = _curAnchor;
    var self = this;
    var _distance = function (anchor, cx, cy, xy, wh, r, tr) {
      var ax = xy[0] + anchor.x * wh[0];
      var ay = xy[1] + anchor.y * wh[1];
      var acx = xy[0] + wh[0] / 2;
      var acy = xy[1] + wh[1] / 2;
      if (r != null && r !== 0) {
        var rotated = jsPlumbUtil.rotatePoint([ax, ay], [acx, acy], r);
        ax = rotated[0];
        ay = rotated[1];
      }
      return (
        Math.sqrt(Math.pow(cx - ax, 2) + Math.pow(cy - ay, 2)) +
        Math.sqrt(Math.pow(acx - ax, 2) + Math.pow(acy - ay, 2))
      );
    };
    var _anchorSelector =
      params$jscomp$0.selector ||
      function (xy, wh, txy, twh, r, tr, anchors) {
        var cx = txy[0] + twh[0] / 2;
        var cy = txy[1] + twh[1] / 2;
        var minIdx = -1;
        var minDist = Infinity;
        for (var i = 0; i < anchors.length; i++) {
          var d = _distance(anchors[i], cx, cy, xy, wh, r, tr);
          if (d < minDist) {
            minIdx = i + 0;
            minDist = d;
          }
        }
        return anchors[minIdx];
      };
    this.compute = function (params) {
      var xy = params.xy;
      var wh = params.wh;
      var txy = params.txy;
      var twh = params.twh;
      var r = params.rotation;
      var tr = params.tRotation;
      this.timestamp = params.timestamp;
      if (this.locked || txy == null || twh == null) {
        this.lastReturnValue = _curAnchor.compute(params);
        return this.lastReturnValue;
      } else params.timestamp = null;
      _curAnchor = _anchorSelector(xy, wh, txy, twh, r, tr, this.anchors);
      this.x = _curAnchor.x;
      this.y = _curAnchor.y;
      if (_curAnchor !== _lastAnchor) this.fire("anchorChanged", _curAnchor);
      _lastAnchor = _curAnchor;
      this.lastReturnValue = _curAnchor.compute(params);
      return this.lastReturnValue;
    };
    this.getCurrentLocation = function (params) {
      return _curAnchor != null ? _curAnchor.getCurrentLocation(params) : null;
    };
    this.getOrientation = function (_endpoint) {
      return _curAnchor != null ? _curAnchor.getOrientation(_endpoint) : [0, 0];
    };
    this.over = function (anchor, endpoint) {
      if (_curAnchor != null) _curAnchor.over(anchor, endpoint);
    };
    this.out = function () {
      if (_curAnchor != null) _curAnchor.out();
    };
    this.setAnchor = function (a) {
      _curAnchor = a;
    };
    this.getCssClass = function () {
      return (_curAnchor && _curAnchor.getCssClass()) || "";
    };
    this.setAnchorCoordinates = function (coords) {
      var idx = jsPlumbUtil.findWithFunction(this.anchors, function (a) {
        return a.x === coords[0] && a.y === coords[1];
      });
      if (idx !== -1) {
        this.setAnchor(this.anchors[idx]);
        return true;
      } else return false;
    };
  };
  _ju.extend(_jp.DynamicAnchor, _jp.Anchor);
  var _curryAnchor = function (x, y, ox, oy, type, fnInit) {
    _jp.Anchors[type] = function (params) {
      var a = params.jsPlumbInstance.makeAnchor(
        [x, y, ox, oy, 0, 0],
        params.elementId,
        params.jsPlumbInstance
      );
      a.type = type;
      if (fnInit) fnInit(a, params);
      return a;
    };
  };
  _curryAnchor(0.5, 0, 0, -1, "TopCenter");
  _curryAnchor(0.5, 1, 0, 1, "BottomCenter");
  _curryAnchor(0, 0.5, -1, 0, "LeftMiddle");
  _curryAnchor(1, 0.5, 1, 0, "RightMiddle");
  _curryAnchor(0.5, 0, 0, -1, "Top");
  _curryAnchor(0.5, 1, 0, 1, "Bottom");
  _curryAnchor(0, 0.5, -1, 0, "Left");
  _curryAnchor(1, 0.5, 1, 0, "Right");
  _curryAnchor(0.5, 0.5, 0, 0, "Center");
  _curryAnchor(1, 0, 0, -1, "TopRight");
  _curryAnchor(1, 1, 0, 1, "BottomRight");
  _curryAnchor(0, 0, 0, -1, "TopLeft");
  _curryAnchor(0, 1, 0, 1, "BottomLeft");
  _jp.Defaults.DynamicAnchors = function (params) {
    return params.jsPlumbInstance.makeAnchors(
      ["TopCenter", "RightMiddle", "BottomCenter", "LeftMiddle"],
      params.elementId,
      params.jsPlumbInstance
    );
  };
  _jp.Anchors.AutoDefault = function (params) {
    var a = params.jsPlumbInstance.makeDynamicAnchor(
      _jp.Defaults.DynamicAnchors(params)
    );
    a.type = "AutoDefault";
    return a;
  };
  var _curryContinuousAnchor = function (type, faces) {
    _jp.Anchors[type] = function (params) {
      var a = params.jsPlumbInstance.makeAnchor(
        ["Continuous", { faces }],
        params.elementId,
        params.jsPlumbInstance
      );
      a.type = type;
      return a;
    };
  };
  _jp.Anchors.Continuous = function (params) {
    return params.jsPlumbInstance.continuousAnchorFactory.get(params);
  };
  _curryContinuousAnchor("ContinuousLeft", ["left"]);
  _curryContinuousAnchor("ContinuousTop", ["top"]);
  _curryContinuousAnchor("ContinuousBottom", ["bottom"]);
  _curryContinuousAnchor("ContinuousRight", ["right"]);
  _curryAnchor(0, 0, 0, 0, "Assign", function (anchor, params) {
    var pf = params.position || "Fixed";
    anchor.positionFinder =
      pf.constructor === String
        ? params.jsPlumbInstance.AnchorPositionFinders[pf]
        : pf;
    anchor.constructorParams = params;
  });
  root.jsPlumbInstance.prototype.AnchorPositionFinders = {
    Fixed: function (dp, ep, es) {
      return [(dp.left - ep.left) / es[0], (dp.top - ep.top) / es[1]];
    },
    Grid: function (dp, ep, es, params) {
      var dx = dp.left - ep.left;
      var dy = dp.top - ep.top;
      var gx = es[0] / params.grid[0];
      var gy = es[1] / params.grid[1];
      var mx = Math.floor(dx / gx);
      var my = Math.floor(dy / gy);
      return [(mx * gx + gx / 2) / es[0], (my * gy + gy / 2) / es[1]];
    },
  };
  _jp.Anchors.Perimeter = function (params$jscomp$0) {
    params$jscomp$0 = params$jscomp$0 || {};
    var anchorCount = params$jscomp$0.anchorCount || 60;
    var shape = params$jscomp$0.shape;
    if (!shape) throw new Error("no shape supplied to Perimeter Anchor type");
    var _circle = function () {
      var r = 0.5;
      var step = (Math.PI * 2) / anchorCount;
      var current = 0;
      var a = [];
      for (var i = 0; i < anchorCount; i++) {
        var x = r + r * Math.sin(current);
        var y = r + r * Math.cos(current);
        a.push([x, y, 0, 0]);
        current += step;
      }
      return a;
    };
    var _path = function (segments) {
      var anchorsPerFace = anchorCount / segments.length;
      var a = [];
      var _computeFace = function (x1, y1, x2, y2, fractionalLength, ox, oy) {
        anchorsPerFace = anchorCount * fractionalLength;
        var dx = (x2 - x1) / anchorsPerFace;
        var dy = (y2 - y1) / anchorsPerFace;
        for (var i = 0; i < anchorsPerFace; i++)
          a.push([
            x1 + dx * i,
            y1 + dy * i,
            ox == null ? 0 : ox,
            oy == null ? 0 : oy,
          ]);
      };
      for (var i$jscomp$0 = 0; i$jscomp$0 < segments.length; i$jscomp$0++)
        _computeFace.apply(null, segments[i$jscomp$0]);
      return a;
    };
    var _shape = function (faces) {
      var s = [];
      for (var i = 0; i < faces.length; i++)
        s.push([
          faces[i][0],
          faces[i][1],
          faces[i][2],
          faces[i][3],
          1 / faces.length,
          faces[i][4],
          faces[i][5],
        ]);
      return _path(s);
    };
    var _rectangle = function () {
      return _shape([
        [0, 0, 1, 0, 0, -1],
        [1, 0, 1, 1, 1, 0],
        [1, 1, 0, 1, 0, 1],
        [0, 1, 0, 0, -1, 0],
      ]);
    };
    var _shapes = {
      Circle: _circle,
      Ellipse: _circle,
      Diamond: function () {
        return _shape([
          [0.5, 0, 1, 0.5],
          [1, 0.5, 0.5, 1],
          [0.5, 1, 0, 0.5],
          [0, 0.5, 0.5, 0],
        ]);
      },
      Rectangle: _rectangle,
      Square: _rectangle,
      Triangle: function () {
        return _shape([
          [0.5, 0, 1, 1],
          [1, 1, 0, 1],
          [0, 1, 0.5, 0],
        ]);
      },
      Path: function (params) {
        var points = params.points;
        var p = [];
        var tl = 0;
        for (var i = 0; i < points.length - 1; i++) {
          var l = Math.sqrt(
            Math.pow(points[i][2] - points[i][0]) +
              Math.pow(points[i][3] - points[i][1])
          );
          tl += l;
          p.push([
            points[i][0],
            points[i][1],
            points[i + 1][0],
            points[i + 1][1],
            l,
          ]);
        }
        for (var j = 0; j < p.length; j++) p[j][4] = p[j][4] / tl;
        return _path(p);
      },
    };
    var _rotate = function (points, amountInDegrees) {
      var o = [];
      var theta = (amountInDegrees / 180) * Math.PI;
      for (var i = 0; i < points.length; i++) {
        var _x = points[i][0] - 0.5;
        var _y = points[i][1] - 0.5;
        o.push([
          0.5 + (_x * Math.cos(theta) - _y * Math.sin(theta)),
          0.5 + (_x * Math.sin(theta) + _y * Math.cos(theta)),
          points[i][2],
          points[i][3],
        ]);
      }
      return o;
    };
    if (!_shapes[shape])
      throw new Error(
        "Shape [" + shape + "] is unknown by Perimeter Anchor type"
      );
    var da = _shapes[shape](params$jscomp$0);
    if (params$jscomp$0.rotation) da = _rotate(da, params$jscomp$0.rotation);
    var a$jscomp$0 = params$jscomp$0.jsPlumbInstance.makeDynamicAnchor(da);
    a$jscomp$0.type = "Perimeter";
    return a$jscomp$0;
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _ju = root.jsPlumbUtil;
  var _jp = root.jsPlumb;
  _jp.DefaultRouter = function (jsPlumbInstance) {
    this.jsPlumbInstance = jsPlumbInstance;
    this.anchorManager = new _jp.AnchorManager({ jsPlumbInstance });
    this.sourceOrTargetChanged = function (
      originalId,
      newId,
      connection,
      newElement,
      anchorIndex
    ) {
      this.anchorManager.sourceOrTargetChanged(
        originalId,
        newId,
        connection,
        newElement,
        anchorIndex
      );
    };
    this.reset = function () {
      this.anchorManager.reset();
    };
    this.changeId = function (oldId, newId) {
      this.anchorManager.changeId(oldId, newId);
    };
    this.elementRemoved = function (elementId) {
      this.anchorManager.elementRemoved(elementId);
    };
    this.newConnection = function (conn) {
      this.anchorManager.newConnection(conn);
    };
    this.connectionDetached = function (connInfo, doNotRedraw) {
      this.anchorManager.connectionDetached(connInfo, doNotRedraw);
    };
    this.redraw = function (
      elementId,
      ui,
      timestamp,
      offsetToUI,
      clearEdits,
      doNotRecalcEndpoint
    ) {
      return this.anchorManager.redraw(
        elementId,
        ui,
        timestamp,
        offsetToUI,
        clearEdits,
        doNotRecalcEndpoint
      );
    };
    this.deleteEndpoint = function (endpoint) {
      this.anchorManager.deleteEndpoint(endpoint);
    };
    this.rehomeEndpoint = function (ep, currentId, element) {
      this.anchorManager.rehomeEndpoint(ep, currentId, element);
    };
    this.addEndpoint = function (endpoint, elementId) {
      this.anchorManager.addEndpoint(endpoint, elementId);
    };
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var _jg = root.Biltong;
  _jp.Segments = {
    AbstractSegment: function (params) {
      this.params = params;
      this.findClosestPointOnPath = function (x, y) {
        return { d: Infinity, x: null, y: null, l: null };
      };
      this.getBounds = function () {
        return {
          minX: Math.min(params.x1, params.x2),
          minY: Math.min(params.y1, params.y2),
          maxX: Math.max(params.x1, params.x2),
          maxY: Math.max(params.y1, params.y2),
        };
      };
      this.lineIntersection = function (x1, y1, x2, y2) {
        return [];
      };
      this.boxIntersection = function (x, y, w, h) {
        var a = [];
        a.push.apply(a, this.lineIntersection(x, y, x + w, y));
        a.push.apply(a, this.lineIntersection(x + w, y, x + w, y + h));
        a.push.apply(a, this.lineIntersection(x + w, y + h, x, y + h));
        a.push.apply(a, this.lineIntersection(x, y + h, x, y));
        return a;
      };
      this.boundingBoxIntersection = function (box) {
        return this.boxIntersection(box.x, box.y, box.w, box.y);
      };
    },
    Straight: function (params) {
      var _super = _jp.Segments.AbstractSegment.apply(this, arguments);
      var length;
      var m;
      var m2;
      var x1;
      var x2;
      var y1;
      var y2;
      var _recalc = function () {
        length = Math.sqrt(Math.pow(x2 - x1, 2) + Math.pow(y2 - y1, 2));
        m = _jg.gradient({ x: x1, y: y1 }, { x: x2, y: y2 });
        m2 = -1 / m;
      };
      this.type = "Straight";
      this.getLength = function () {
        return length;
      };
      this.getGradient = function () {
        return m;
      };
      this.getCoordinates = function () {
        return { x1, y1, x2, y2 };
      };
      this.setCoordinates = function (coords) {
        x1 = coords.x1;
        y1 = coords.y1;
        x2 = coords.x2;
        y2 = coords.y2;
        _recalc();
      };
      this.setCoordinates({
        x1: params.x1,
        y1: params.y1,
        x2: params.x2,
        y2: params.y2,
      });
      this.getBounds = function () {
        return {
          minX: Math.min(x1, x2),
          minY: Math.min(y1, y2),
          maxX: Math.max(x1, x2),
          maxY: Math.max(y1, y2),
        };
      };
      this.pointOnPath = function (location, absolute) {
        if (location === 0 && !absolute) return { x: x1, y: y1 };
        else if (location === 1 && !absolute) return { x: x2, y: y2 };
        else {
          var l = absolute
            ? location > 0
              ? location
              : length + location
            : location * length;
          return _jg.pointOnLine({ x: x1, y: y1 }, { x: x2, y: y2 }, l);
        }
      };
      this.gradientAtPoint = function (_) {
        return m;
      };
      this.pointAlongPathFrom = function (location, distance, absolute) {
        var p = this.pointOnPath(location, absolute);
        var farAwayPoint = distance <= 0 ? { x: x1, y: y1 } : { x: x2, y: y2 };
        if (distance <= 0 && Math.abs(distance) > 1) distance *= -1;
        return _jg.pointOnLine(p, farAwayPoint, distance);
      };
      var within = function (a, b, c) {
        return c >= Math.min(a, b) && c <= Math.max(a, b);
      };
      var closest = function (a, b, c) {
        return Math.abs(c - a) < Math.abs(c - b) ? a : b;
      };
      this.findClosestPointOnPath = function (x, y) {
        var out = { d: Infinity, x: null, y: null, l: null, x1, x2, y1, y2 };
        if (m === 0) {
          out.y = y1;
          out.x = within(x1, x2, x) ? x : closest(x1, x2, x);
        } else if (m === Infinity || m === -Infinity) {
          out.x = x1;
          out.y = within(y1, y2, y) ? y : closest(y1, y2, y);
        } else {
          var b = y1 - m * x1;
          var b2 = y - m2 * x;
          var _x1 = (b2 - b) / (m - m2);
          var _y1 = m * _x1 + b;
          out.x = within(x1, x2, _x1) ? _x1 : closest(x1, x2, _x1);
          out.y = within(y1, y2, _y1) ? _y1 : closest(y1, y2, _y1);
        }
        var fractionInSegment = _jg.lineLength([out.x, out.y], [x1, y1]);
        out.d = _jg.lineLength([x, y], [out.x, out.y]);
        out.l = fractionInSegment / length;
        return out;
      };
      var _pointLiesBetween = function (q, p1, p2) {
        return p2 > p1 ? p1 <= q && q <= p2 : p1 >= q && q >= p2;
      };
      var _plb = _pointLiesBetween;
      this.lineIntersection = function (_x1, _y1, _x2, _y2) {
        var m2 = Math.abs(_jg.gradient({ x: _x1, y: _y1 }, { x: _x2, y: _y2 }));
        var m1 = Math.abs(m);
        var b = m1 === Infinity ? x1 : y1 - m1 * x1;
        var out = [];
        var b2 = m2 === Infinity ? _x1 : _y1 - m2 * _x1;
        if (m2 !== m1)
          if (m2 === Infinity && m1 === 0) {
            if (_plb(_x1, x1, x2) && _plb(y1, _y1, _y2)) out = [_x1, y1];
          } else if (m2 === 0 && m1 === Infinity) {
            if (_plb(_y1, y1, y2) && _plb(x1, _x1, _x2)) out = [x1, _y1];
          } else if (m2 === Infinity) {
            var X = _x1;
            if (_plb(X, x1, x2)) {
              var Y = m1 * _x1 + b;
              if (_plb(Y, _y1, _y2)) out = [X, Y];
            }
          } else if (m2 === 0) {
            Y = _y1;
            if (_plb(Y, y1, y2)) {
              X = (_y1 - b) / m1;
              if (_plb(X, _x1, _x2)) out = [X, Y];
            }
          } else {
            X = (b2 - b) / (m1 - m2);
            Y = m1 * X + b;
            if (_plb(X, x1, x2) && _plb(Y, y1, y2)) out = [X, Y];
          }
        return out;
      };
      this.boxIntersection = function (x, y, w, h) {
        var a = [];
        a.push.apply(a, this.lineIntersection(x, y, x + w, y));
        a.push.apply(a, this.lineIntersection(x + w, y, x + w, y + h));
        a.push.apply(a, this.lineIntersection(x + w, y + h, x, y + h));
        a.push.apply(a, this.lineIntersection(x, y + h, x, y));
        return a;
      };
      this.boundingBoxIntersection = function (box) {
        return this.boxIntersection(box.x, box.y, box.w, box.h);
      };
    },
    Arc: function (params) {
      var _super = _jp.Segments.AbstractSegment.apply(this, arguments);
      var _calcAngle = function (_x, _y) {
        return _jg.theta([params.cx, params.cy], [_x, _y]);
      };
      var _calcAngleForLocation = function (segment, location) {
        if (segment.anticlockwise) {
          var sa =
            segment.startAngle < segment.endAngle
              ? segment.startAngle + TWO_PI
              : segment.startAngle;
          var s = Math.abs(sa - segment.endAngle);
          return sa - s * location;
        } else {
          var ea =
            segment.endAngle < segment.startAngle
              ? segment.endAngle + TWO_PI
              : segment.endAngle;
          var ss = Math.abs(ea - segment.startAngle);
          return segment.startAngle + ss * location;
        }
      };
      var TWO_PI = 2 * Math.PI;
      this.radius = params.r;
      this.anticlockwise = params.ac;
      this.type = "Arc";
      if (params.startAngle && params.endAngle) {
        this.startAngle = params.startAngle;
        this.endAngle = params.endAngle;
        this.x1 = params.cx + this.radius * Math.cos(params.startAngle);
        this.y1 = params.cy + this.radius * Math.sin(params.startAngle);
        this.x2 = params.cx + this.radius * Math.cos(params.endAngle);
        this.y2 = params.cy + this.radius * Math.sin(params.endAngle);
      } else {
        this.startAngle = _calcAngle(params.x1, params.y1);
        this.endAngle = _calcAngle(params.x2, params.y2);
        this.x1 = params.x1;
        this.y1 = params.y1;
        this.x2 = params.x2;
        this.y2 = params.y2;
      }
      if (this.endAngle < 0) this.endAngle += TWO_PI;
      if (this.startAngle < 0) this.startAngle += TWO_PI;
      var ea =
        this.endAngle < this.startAngle
          ? this.endAngle + TWO_PI
          : this.endAngle;
      this.sweep = Math.abs(ea - this.startAngle);
      if (this.anticlockwise) this.sweep = TWO_PI - this.sweep;
      var circumference = 2 * Math.PI * this.radius;
      var frac = this.sweep / TWO_PI;
      var length = circumference * frac;
      this.getLength = function () {
        return length;
      };
      this.getBounds = function () {
        return {
          minX: params.cx - params.r,
          maxX: params.cx + params.r,
          minY: params.cy - params.r,
          maxY: params.cy + params.r,
        };
      };
      var VERY_SMALL_VALUE = 1e-10;
      var gentleRound = function (n) {
        var f = Math.floor(n);
        var r = Math.ceil(n);
        if (n - f < VERY_SMALL_VALUE) return f;
        else if (r - n < VERY_SMALL_VALUE) return r;
        return n;
      };
      this.pointOnPath = function (location, absolute) {
        if (location === 0)
          return { x: this.x1, y: this.y1, theta: this.startAngle };
        else if (location === 1)
          return { x: this.x2, y: this.y2, theta: this.endAngle };
        if (absolute) location /= length;
        var angle = _calcAngleForLocation(this, location);
        var _x = params.cx + params.r * Math.cos(angle);
        var _y = params.cy + params.r * Math.sin(angle);
        return { x: gentleRound(_x), y: gentleRound(_y), theta: angle };
      };
      this.gradientAtPoint = function (location, absolute) {
        var p = this.pointOnPath(location, absolute);
        var m = _jg.normal([params.cx, params.cy], [p.x, p.y]);
        if (!this.anticlockwise && (m === Infinity || m === -Infinity)) m *= -1;
        return m;
      };
      this.pointAlongPathFrom = function (location, distance, absolute) {
        var p = this.pointOnPath(location, absolute);
        var arcSpan = (distance / circumference) * 2 * Math.PI;
        var dir = this.anticlockwise ? -1 : 1;
        var startAngle = p.theta + dir * arcSpan;
        var startX = params.cx + this.radius * Math.cos(startAngle);
        var startY = params.cy + this.radius * Math.sin(startAngle);
        return { x: startX, y: startY };
      };
    },
    Bezier: function (params) {
      this.curve = [
        { x: params.x1, y: params.y1 },
        { x: params.cp1x, y: params.cp1y },
        { x: params.cp2x, y: params.cp2y },
        { x: params.x2, y: params.y2 },
      ];
      var _isPoint = function (c) {
        return c[0].x === c[1].x && c[0].y === c[1].y;
      };
      var _dist = function (p1, p2) {
        return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
      };
      var _compute = function (loc) {
        var EMPTY_POINT = { x: 0, y: 0 };
        if (loc === 0) return this.curve[0];
        var degree = this.curve.length - 1;
        if (loc === 1) return this.curve[degree];
        var o = this.curve;
        var s = 1 - loc;
        if (degree === 0) return this.curve[0];
        if (degree === 1)
          return { x: s * o[0].x + loc * o[1].x, y: s * o[0].y + loc * o[1].y };
        if (degree < 4) {
          var l = s * s;
          var h = loc * loc;
          var u = 0;
          if (degree === 2) {
            o = [o[0], o[1], o[2], EMPTY_POINT];
            var m = l;
            var g = 2 * (s * loc);
            var f = h;
          } else if (degree === 3) {
            m = l * s;
            g = 3 * (l * loc);
            f = 3 * (s * h);
            u = loc * h;
          }
          return {
            x: m * o[0].x + g * o[1].x + f * o[2].x + u * o[3].x,
            y: m * o[0].y + g * o[1].y + f * o[2].y + u * o[3].y,
          };
        } else return EMPTY_POINT;
      }.bind(this);
      var _getLUT = function (steps) {
        var out = [];
        steps--;
        for (var n = 0; n <= steps; n++) out.push(_compute(n / steps));
        return out;
      };
      var _computeLength = function () {
        if (_isPoint(this.curve)) this.length = 0;
        var steps = 16;
        var lut = _getLUT(steps);
        this.length = 0;
        for (var i = 0; i < steps - 1; i++) {
          var a = lut[i];
          var b = lut[i + 1];
          this.length += _dist(a, b);
        }
      }.bind(this);
      var _super = _jp.Segments.AbstractSegment.apply(this, arguments);
      this.bounds = {
        minX: Math.min(params.x1, params.x2, params.cp1x, params.cp2x),
        minY: Math.min(params.y1, params.y2, params.cp1y, params.cp2y),
        maxX: Math.max(params.x1, params.x2, params.cp1x, params.cp2x),
        maxY: Math.max(params.y1, params.y2, params.cp1y, params.cp2y),
      };
      this.type = "Bezier";
      _computeLength();
      var _translateLocation = function (_curve, location, absolute) {
        if (absolute)
          location = root.jsBezier.locationAlongCurveFrom(
            _curve,
            location > 0 ? 0 : 1,
            location
          );
        return location;
      };
      this.pointOnPath = function (location, absolute) {
        location = _translateLocation(this.curve, location, absolute);
        return root.jsBezier.pointOnCurve(this.curve, location);
      };
      this.gradientAtPoint = function (location, absolute) {
        location = _translateLocation(this.curve, location, absolute);
        return root.jsBezier.gradientAtPoint(this.curve, location);
      };
      this.pointAlongPathFrom = function (location, distance, absolute) {
        location = _translateLocation(this.curve, location, absolute);
        return root.jsBezier.pointAlongCurveFrom(
          this.curve,
          location,
          distance
        );
      };
      this.getLength = function () {
        return this.length;
      };
      this.getBounds = function () {
        return this.bounds;
      };
      this.findClosestPointOnPath = function (x, y) {
        var p = root.jsBezier.nearestPointOnCurve({ x, y }, this.curve);
        return {
          d: Math.sqrt(Math.pow(p.point.x - x, 2) + Math.pow(p.point.y - y, 2)),
          x: p.point.x,
          y: p.point.y,
          l: 1 - p.location,
          s: this,
        };
      };
      this.lineIntersection = function (x1, y1, x2, y2) {
        return root.jsBezier.lineIntersection(x1, y1, x2, y2, this.curve);
      };
    },
  };
  _jp.SegmentRenderer = {
    getPath: function (segment, isFirstSegment) {
      return {
        Straight: function (isFirstSegment) {
          var d = segment.getCoordinates();
          return (
            (isFirstSegment ? "M " + d.x1 + " " + d.y1 + " " : "") +
            "L " +
            d.x2 +
            " " +
            d.y2
          );
        },
        Bezier: function (isFirstSegment) {
          var d = segment.params;
          return (
            (isFirstSegment ? "M " + d.x2 + " " + d.y2 + " " : "") +
            "C " +
            d.cp2x +
            " " +
            d.cp2y +
            " " +
            d.cp1x +
            " " +
            d.cp1y +
            " " +
            d.x1 +
            " " +
            d.y1
          );
        },
        Arc: function (isFirstSegment) {
          var d = segment.params;
          var laf = segment.sweep > Math.PI ? 1 : 0;
          var sf = segment.anticlockwise ? 0 : 1;
          return (
            (isFirstSegment ? "M" + segment.x1 + " " + segment.y1 + " " : "") +
            "A " +
            segment.radius +
            " " +
            d.r +
            " 0 " +
            laf +
            "," +
            sf +
            " " +
            segment.x2 +
            " " +
            segment.y2
          );
        },
      }[segment.type](isFirstSegment);
    },
  };
  var AbstractComponent = function () {
    this.resetBounds = function () {
      this.bounds = {
        minX: Infinity,
        minY: Infinity,
        maxX: -Infinity,
        maxY: -Infinity,
      };
    };
    this.resetBounds();
  };
  _jp.Connectors.AbstractConnector = function (params$jscomp$0) {
    AbstractComponent.apply(this, arguments);
    var segments = [];
    var totalLength = 0;
    var segmentProportions = [];
    var segmentProportionalLengths = [];
    var stub = params$jscomp$0.stub || 0;
    var sourceStub = _ju.isArray(stub) ? stub[0] : stub;
    var targetStub = _ju.isArray(stub) ? stub[1] : stub;
    var gap = params$jscomp$0.gap || 0;
    var sourceGap = _ju.isArray(gap) ? gap[0] : gap;
    var targetGap = _ju.isArray(gap) ? gap[1] : gap;
    var userProvidedSegments = null;
    var paintInfo = null;
    this.getPathData = function () {
      var p = "";
      for (var i = 0; i < segments.length; i++) {
        p += _jp.SegmentRenderer.getPath(segments[i], i === 0);
        p += " ";
      }
      return p;
    };
    this.findSegmentForPoint = function (x, y) {
      var out = { d: Infinity, s: null, x: null, y: null, l: null };
      for (var i = 0; i < segments.length; i++) {
        var _s = segments[i].findClosestPointOnPath(x, y);
        if (_s.d < out.d) {
          out.d = _s.d;
          out.l = _s.l;
          out.x = _s.x;
          out.y = _s.y;
          out.s = segments[i];
          out.x1 = _s.x1;
          out.x2 = _s.x2;
          out.y1 = _s.y1;
          out.y2 = _s.y2;
          out.index = i;
          out.connectorLocation =
            segmentProportions[i][0] +
            _s.l * (segmentProportions[i][1] - segmentProportions[i][0]);
        }
      }
      return out;
    };
    this.lineIntersection = function (x1, y1, x2, y2) {
      var out = [];
      for (var i = 0; i < segments.length; i++)
        out.push.apply(out, segments[i].lineIntersection(x1, y1, x2, y2));
      return out;
    };
    this.boxIntersection = function (x, y, w, h) {
      var out = [];
      for (var i = 0; i < segments.length; i++)
        out.push.apply(out, segments[i].boxIntersection(x, y, w, h));
      return out;
    };
    this.boundingBoxIntersection = function (box) {
      var out = [];
      for (var i = 0; i < segments.length; i++)
        out.push.apply(out, segments[i].boundingBoxIntersection(box));
      return out;
    };
    var _updateSegmentProportions = function () {
      var curLoc = 0;
      for (var i = 0; i < segments.length; i++) {
        var sl = segments[i].getLength();
        segmentProportionalLengths[i] = sl / totalLength;
        segmentProportions[i] = [curLoc, (curLoc += sl / totalLength)];
      }
    };
    var _findSegmentForLocation = function (location, absolute) {
      var i;
      if (absolute)
        location =
          location > 0
            ? location / totalLength
            : (totalLength + location) / totalLength;
      if (location === 1) {
        var idx = segments.length - 1;
        var inSegmentProportion = 1;
      } else if (location === 0) {
        inSegmentProportion = 0;
        idx = 0;
      } else if (location >= 0.5) {
        idx = 0;
        inSegmentProportion = 0;
        for (i = segmentProportions.length - 1; i > -1; i--)
          if (
            segmentProportions[i][1] >= location &&
            segmentProportions[i][0] <= location
          ) {
            idx = i;
            inSegmentProportion =
              (location - segmentProportions[i][0]) /
              segmentProportionalLengths[i];
            break;
          }
      } else {
        idx = segmentProportions.length - 1;
        inSegmentProportion = 1;
        for (i = 0; i < segmentProportions.length; i++)
          if (segmentProportions[i][1] >= location) {
            idx = i;
            inSegmentProportion =
              (location - segmentProportions[i][0]) /
              segmentProportionalLengths[i];
            break;
          }
      }
      return {
        segment: segments[idx],
        proportion: inSegmentProportion,
        index: idx,
      };
    };
    var _addSegment = function (conn, type, params) {
      if (params.x1 === params.x2 && params.y1 === params.y2) return;
      var s = new _jp.Segments[type](params);
      segments.push(s);
      totalLength += s.getLength();
      conn.updateBounds(s);
    };
    var _clearSegments = function () {
      totalLength =
        segments.length =
        segmentProportions.length =
        segmentProportionalLengths.length =
          0;
    };
    this.setSegments = function (_segs) {
      userProvidedSegments = [];
      totalLength = 0;
      for (var i = 0; i < _segs.length; i++) {
        userProvidedSegments.push(_segs[i]);
        totalLength += _segs[i].getLength();
      }
    };
    this.getLength = function () {
      return totalLength;
    };
    var _prepareCompute = function (params) {
      this.strokeWidth = params.strokeWidth;
      var segment = _jg.quadrant(params.sourcePos, params.targetPos);
      var swapX = params.targetPos[0] < params.sourcePos[0];
      var swapY = params.targetPos[1] < params.sourcePos[1];
      var lw = params.strokeWidth || 1;
      var so = params.sourceEndpoint.anchor.getOrientation(
        params.sourceEndpoint
      );
      var to = params.targetEndpoint.anchor.getOrientation(
        params.targetEndpoint
      );
      var x = swapX ? params.targetPos[0] : params.sourcePos[0];
      var y = swapY ? params.targetPos[1] : params.sourcePos[1];
      var w = Math.abs(params.targetPos[0] - params.sourcePos[0]);
      var h = Math.abs(params.targetPos[1] - params.sourcePos[1]);
      if ((so[0] === 0 && so[1] === 0) || (to[0] === 0 && to[1] === 0)) {
        var index = w > h ? 0 : 1;
        var oIndex = [1, 0][index];
        so = [];
        to = [];
        so[index] = params.sourcePos[index] > params.targetPos[index] ? -1 : 1;
        to[index] = params.sourcePos[index] > params.targetPos[index] ? 1 : -1;
        so[oIndex] = 0;
        to[oIndex] = 0;
      }
      var sx = swapX ? w + sourceGap * so[0] : sourceGap * so[0];
      var sy = swapY ? h + sourceGap * so[1] : sourceGap * so[1];
      var tx = swapX ? targetGap * to[0] : w + targetGap * to[0];
      var ty = swapY ? targetGap * to[1] : h + targetGap * to[1];
      var oProduct = so[0] * to[0] + so[1] * to[1];
      var result = {
        sx,
        sy,
        tx,
        ty,
        lw,
        xSpan: Math.abs(tx - sx),
        ySpan: Math.abs(ty - sy),
        mx: (sx + tx) / 2,
        my: (sy + ty) / 2,
        so,
        to,
        x,
        y,
        w,
        h,
        segment,
        startStubX: sx + so[0] * sourceStub,
        startStubY: sy + so[1] * sourceStub,
        endStubX: tx + to[0] * targetStub,
        endStubY: ty + to[1] * targetStub,
        isXGreaterThanStubTimes2: Math.abs(sx - tx) > sourceStub + targetStub,
        isYGreaterThanStubTimes2: Math.abs(sy - ty) > sourceStub + targetStub,
        opposite: oProduct === -1,
        perpendicular: oProduct === 0,
        orthogonal: oProduct === 1,
        sourceAxis: so[0] === 0 ? "y" : "x",
        points: [x, y, w, h, sx, sy, tx, ty],
        stubs: [sourceStub, targetStub],
      };
      result.anchorOrientation = result.opposite
        ? "opposite"
        : result.orthogonal
        ? "orthogonal"
        : "perpendicular";
      return result;
    };
    this.getSegments = function () {
      return segments;
    };
    this.updateBounds = function (segment) {
      var segBounds = segment.getBounds();
      this.bounds.minX = Math.min(this.bounds.minX, segBounds.minX);
      this.bounds.maxX = Math.max(this.bounds.maxX, segBounds.maxX);
      this.bounds.minY = Math.min(this.bounds.minY, segBounds.minY);
      this.bounds.maxY = Math.max(this.bounds.maxY, segBounds.maxY);
    };
    var dumpSegmentsToConsole = function () {
      console.log("SEGMENTS:");
      for (var i = 0; i < segments.length; i++)
        console.log(
          segments[i].type,
          segments[i].getLength(),
          segmentProportions[i]
        );
    };
    this.pointOnPath = function (location, absolute) {
      var seg = _findSegmentForLocation(location, absolute);
      return (
        (seg.segment && seg.segment.pointOnPath(seg.proportion, false)) || [
          0, 0,
        ]
      );
    };
    this.gradientAtPoint = function (location, absolute) {
      var seg = _findSegmentForLocation(location, absolute);
      return (
        (seg.segment && seg.segment.gradientAtPoint(seg.proportion, false)) || 0
      );
    };
    this.pointAlongPathFrom = function (location, distance, absolute) {
      var seg = _findSegmentForLocation(location, absolute);
      return (
        (seg.segment &&
          seg.segment.pointAlongPathFrom(seg.proportion, distance, false)) || [
          0, 0,
        ]
      );
    };
    this.compute = function (params) {
      paintInfo = _prepareCompute.call(this, params);
      _clearSegments();
      this._compute(paintInfo, params);
      this.x = paintInfo.points[0];
      this.y = paintInfo.points[1];
      this.w = paintInfo.points[2];
      this.h = paintInfo.points[3];
      this.segment = paintInfo.segment;
      _updateSegmentProportions();
    };
    return {
      addSegment: _addSegment,
      prepareCompute: _prepareCompute,
      sourceStub,
      targetStub,
      maxStub: Math.max(sourceStub, targetStub),
      sourceGap,
      targetGap,
      maxGap: Math.max(sourceGap, targetGap),
    };
  };
  _ju.extend(_jp.Connectors.AbstractConnector, AbstractComponent);
  _jp.Endpoints.AbstractEndpoint = function (params) {
    AbstractComponent.apply(this, arguments);
    var compute = (this.compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      var out = this._compute.apply(this, arguments);
      this.x = out[0];
      this.y = out[1];
      this.w = out[2];
      this.h = out[3];
      this.bounds.minX = this.x;
      this.bounds.minY = this.y;
      this.bounds.maxX = this.x + this.w;
      this.bounds.maxY = this.y + this.h;
      return out;
    });
    return { compute, cssClass: params.cssClass };
  };
  _ju.extend(_jp.Endpoints.AbstractEndpoint, AbstractComponent);
  _jp.Endpoints.Dot = function (params) {
    this.type = "Dot";
    var _super = _jp.Endpoints.AbstractEndpoint.apply(this, arguments);
    params = params || {};
    this.radius = params.radius || 10;
    this.defaultOffset = 0.5 * this.radius;
    this.defaultInnerRadius = this.radius / 3;
    this._compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      this.radius = endpointStyle.radius || this.radius;
      var x = anchorPoint[0] - this.radius;
      var y = anchorPoint[1] - this.radius;
      var w = this.radius * 2;
      var h = this.radius * 2;
      if (endpointStyle.stroke) {
        var lw = endpointStyle.strokeWidth || 1;
        x -= lw;
        y -= lw;
        w += lw * 2;
        h += lw * 2;
      }
      return [x, y, w, h, this.radius];
    };
  };
  _ju.extend(_jp.Endpoints.Dot, _jp.Endpoints.AbstractEndpoint);
  _jp.Endpoints.Rectangle = function (params) {
    this.type = "Rectangle";
    var _super = _jp.Endpoints.AbstractEndpoint.apply(this, arguments);
    params = params || {};
    this.width = params.width || 20;
    this.height = params.height || 20;
    this._compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      var width = endpointStyle.width || this.width;
      var height = endpointStyle.height || this.height;
      var x = anchorPoint[0] - width / 2;
      var y = anchorPoint[1] - height / 2;
      return [x, y, width, height];
    };
  };
  _ju.extend(_jp.Endpoints.Rectangle, _jp.Endpoints.AbstractEndpoint);
  var DOMElementEndpoint = function (params) {
    _jp.jsPlumbUIComponent.apply(this, arguments);
    this._jsPlumb.displayElements = [];
  };
  _ju.extend(DOMElementEndpoint, _jp.jsPlumbUIComponent, {
    getDisplayElements: function () {
      return this._jsPlumb.displayElements;
    },
    appendDisplayElement: function (el) {
      this._jsPlumb.displayElements.push(el);
    },
  });
  _jp.Endpoints.Image = function (params) {
    this.type = "Image";
    DOMElementEndpoint.apply(this, arguments);
    _jp.Endpoints.AbstractEndpoint.apply(this, arguments);
    var _onload = params.onload;
    var src = params.src || params.url;
    var clazz = params.cssClass ? " " + params.cssClass : "";
    this._jsPlumb.img = new Image();
    this._jsPlumb.ready = false;
    this._jsPlumb.initialized = false;
    this._jsPlumb.deleted = false;
    this._jsPlumb.widthToUse = params.width;
    this._jsPlumb.heightToUse = params.height;
    this._jsPlumb.endpoint = params.endpoint;
    this._jsPlumb.img.onload = function () {
      if (this._jsPlumb != null) {
        this._jsPlumb.ready = true;
        this._jsPlumb.widthToUse =
          this._jsPlumb.widthToUse || this._jsPlumb.img.width;
        this._jsPlumb.heightToUse =
          this._jsPlumb.heightToUse || this._jsPlumb.img.height;
        if (_onload) _onload(this);
      }
    }.bind(this);
    this._jsPlumb.endpoint.setImage = function (_img, onload) {
      var s = _img.constructor === String ? _img : _img.src;
      _onload = onload;
      this._jsPlumb.img.src = s;
      if (this.canvas != null)
        this.canvas.setAttribute("src", this._jsPlumb.img.src);
    }.bind(this);
    this._jsPlumb.endpoint.setImage(src, _onload);
    this._compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      this.anchorPoint = anchorPoint;
      if (this._jsPlumb.ready)
        return [
          anchorPoint[0] - this._jsPlumb.widthToUse / 2,
          anchorPoint[1] - this._jsPlumb.heightToUse / 2,
          this._jsPlumb.widthToUse,
          this._jsPlumb.heightToUse,
        ];
      else return [0, 0, 0, 0];
    };
    this.canvas = _jp.createElement(
      "img",
      { position: "absolute", margin: 0, padding: 0, outline: 0 },
      this._jsPlumb.instance.endpointClass + clazz
    );
    if (this._jsPlumb.widthToUse)
      this.canvas.setAttribute("width", this._jsPlumb.widthToUse);
    if (this._jsPlumb.heightToUse)
      this.canvas.setAttribute("height", this._jsPlumb.heightToUse);
    this._jsPlumb.instance.appendElement(this.canvas);
    this.actuallyPaint = function (d, style, anchor) {
      if (!this._jsPlumb.deleted) {
        if (!this._jsPlumb.initialized) {
          this.canvas.setAttribute("src", this._jsPlumb.img.src);
          this.appendDisplayElement(this.canvas);
          this._jsPlumb.initialized = true;
        }
        var x = this.anchorPoint[0] - this._jsPlumb.widthToUse / 2;
        var y = this.anchorPoint[1] - this._jsPlumb.heightToUse / 2;
        _ju.sizeElement(
          this.canvas,
          x,
          y,
          this._jsPlumb.widthToUse,
          this._jsPlumb.heightToUse
        );
      }
    };
    this.paint = function (style, anchor) {
      if (this._jsPlumb != null)
        if (this._jsPlumb.ready) this.actuallyPaint(style, anchor);
        else
          root.setTimeout(
            function () {
              this.paint(style, anchor);
            }.bind(this),
            200
          );
    };
  };
  _ju.extend(
    _jp.Endpoints.Image,
    [DOMElementEndpoint, _jp.Endpoints.AbstractEndpoint],
    {
      cleanup: function (force) {
        if (force) {
          this._jsPlumb.deleted = true;
          if (this.canvas) this.canvas.parentNode.removeChild(this.canvas);
          this.canvas = null;
        }
      },
    }
  );
  _jp.Endpoints.Blank = function (params) {
    var _super = _jp.Endpoints.AbstractEndpoint.apply(this, arguments);
    this.type = "Blank";
    DOMElementEndpoint.apply(this, arguments);
    this._compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      return [anchorPoint[0], anchorPoint[1], 10, 0];
    };
    var clazz = params.cssClass ? " " + params.cssClass : "";
    this.canvas = _jp.createElement(
      "div",
      {
        display: "block",
        width: "1px",
        height: "1px",
        background: "transparent",
        position: "absolute",
      },
      this._jsPlumb.instance.endpointClass + clazz
    );
    this._jsPlumb.instance.appendElement(this.canvas);
    this.paint = function (style, anchor) {
      _ju.sizeElement(this.canvas, this.x, this.y, this.w, this.h);
    };
  };
  _ju.extend(
    _jp.Endpoints.Blank,
    [_jp.Endpoints.AbstractEndpoint, DOMElementEndpoint],
    {
      cleanup: function () {
        if (this.canvas && this.canvas.parentNode)
          this.canvas.parentNode.removeChild(this.canvas);
      },
    }
  );
  _jp.Endpoints.Triangle = function (params) {
    this.type = "Triangle";
    _jp.Endpoints.AbstractEndpoint.apply(this, arguments);
    var self = this;
    params = params || {};
    params.width = params.width || 55;
    params.height = params.height || 55;
    this.width = params.width;
    this.height = params.height;
    this._compute = function (
      anchorPoint,
      orientation,
      endpointStyle,
      connectorPaintStyle
    ) {
      var width = endpointStyle.width || self.width;
      var height = endpointStyle.height || self.height;
      var x = anchorPoint[0] - width / 2;
      var y = anchorPoint[1] - height / 2;
      return [x, y, width, height];
    };
  };
  var AbstractOverlay = (_jp.Overlays.AbstractOverlay = function (params) {
    this.visible = true;
    this.isAppendedAtTopLevel = true;
    this.component = params.component;
    this.loc = params.location == null ? 0.5 : params.location;
    this.endpointLoc =
      params.endpointLocation == null ? [0.5, 0.5] : params.endpointLocation;
    this.visible = params.visible !== false;
  });
  AbstractOverlay.prototype = {
    cleanup: function (force) {
      if (force) {
        this.component = null;
        this.canvas = null;
        this.endpointLoc = null;
      }
    },
    reattach: function (instance, component) {},
    setVisible: function (val) {
      this.visible = val;
      this.component.repaint();
    },
    isVisible: function () {
      return this.visible;
    },
    hide: function () {
      this.setVisible(false);
    },
    show: function () {
      this.setVisible(true);
    },
    incrementLocation: function (amount) {
      this.loc += amount;
      this.component.repaint();
    },
    setLocation: function (l) {
      this.loc = l;
      this.component.repaint();
    },
    getLocation: function () {
      return this.loc;
    },
    updateFrom: function () {},
  };
  _jp.Overlays.Arrow = function (params) {
    this.type = "Arrow";
    AbstractOverlay.apply(this, arguments);
    this.isAppendedAtTopLevel = false;
    params = params || {};
    var self = this;
    this.length = params.length || 20;
    this.width = params.width || 20;
    this.id = params.id;
    this.direction = (params.direction || 1) < 0 ? -1 : 1;
    var paintStyle = params.paintStyle || { "stroke-width": 1 };
    var foldback = params.foldback || 0.623;
    this.computeMaxSize = function () {
      return self.width * 1.5;
    };
    this.elementCreated = function (p, component) {
      this.path = p;
      if (params.events)
        for (var i in params.events) _jp.on(p, i, params.events[i]);
    };
    this.draw = function (component, currentConnectionPaintStyle) {
      if (component.pointAlongPathFrom) {
        if (_ju.isString(this.loc) || this.loc > 1 || this.loc < 0) {
          var l = parseInt(this.loc, 10);
          var fromLoc = this.loc < 0 ? 1 : 0;
          var hxy = component.pointAlongPathFrom(fromLoc, l, false);
          var mid = component.pointAlongPathFrom(
            fromLoc,
            l - (this.direction * this.length) / 2,
            false
          );
          var txy = _jg.pointOnLine(hxy, mid, this.length);
        } else if (this.loc === 1) {
          hxy = component.pointOnPath(this.loc);
          mid = component.pointAlongPathFrom(this.loc, -this.length);
          txy = _jg.pointOnLine(hxy, mid, this.length);
          if (this.direction === -1) {
            var _ = txy;
            txy = hxy;
            hxy = _;
          }
        } else if (this.loc === 0) {
          txy = component.pointOnPath(this.loc);
          mid = component.pointAlongPathFrom(this.loc, this.length);
          hxy = _jg.pointOnLine(txy, mid, this.length);
          if (this.direction === -1) {
            var __ = txy;
            txy = hxy;
            hxy = __;
          }
        } else {
          hxy = component.pointAlongPathFrom(
            this.loc,
            (this.direction * this.length) / 2
          );
          mid = component.pointOnPath(this.loc);
          txy = _jg.pointOnLine(hxy, mid, this.length);
        }
        var tail = _jg.perpendicularLineTo(hxy, txy, this.width);
        var cxy = _jg.pointOnLine(hxy, txy, foldback * this.length);
        var d = { hxy, tail, cxy };
        var stroke = paintStyle.stroke || currentConnectionPaintStyle.stroke;
        var fill = paintStyle.fill || currentConnectionPaintStyle.stroke;
        var lineWidth =
          paintStyle.strokeWidth || currentConnectionPaintStyle.strokeWidth;
        return {
          component,
          d,
          "stroke-width": lineWidth,
          stroke,
          fill,
          minX: Math.min(hxy.x, tail[0].x, tail[1].x),
          maxX: Math.max(hxy.x, tail[0].x, tail[1].x),
          minY: Math.min(hxy.y, tail[0].y, tail[1].y),
          maxY: Math.max(hxy.y, tail[0].y, tail[1].y),
        };
      } else return { component, minX: 0, maxX: 0, minY: 0, maxY: 0 };
    };
  };
  _ju.extend(_jp.Overlays.Arrow, AbstractOverlay, {
    updateFrom: function (d) {
      this.length = d.length || this.length;
      this.width = d.width || this.width;
      this.direction = d.direction != null ? d.direction : this.direction;
      this.foldback = d.foldback || this.foldback;
    },
    cleanup: function () {
      if (this.path && this.path.parentNode)
        this.path.parentNode.removeChild(this.path);
    },
  });
  _jp.Overlays.PlainArrow = function (params) {
    params = params || {};
    var p = _jp.extend(params, { foldback: 1 });
    _jp.Overlays.Arrow.call(this, p);
    this.type = "PlainArrow";
  };
  _ju.extend(_jp.Overlays.PlainArrow, _jp.Overlays.Arrow);
  _jp.Overlays.Diamond = function (params) {
    params = params || {};
    var l = params.length || 40;
    var p = _jp.extend(params, { length: l / 2, foldback: 2 });
    _jp.Overlays.Arrow.call(this, p);
    this.type = "Diamond";
  };
  _ju.extend(_jp.Overlays.Diamond, _jp.Overlays.Arrow);
  var _getDimensions = function (component, forceRefresh) {
    if (component._jsPlumb.cachedDimensions == null || forceRefresh)
      component._jsPlumb.cachedDimensions = component.getDimensions();
    return component._jsPlumb.cachedDimensions;
  };
  var AbstractDOMOverlay = function (params) {
    _jp.jsPlumbUIComponent.apply(this, arguments);
    AbstractOverlay.apply(this, arguments);
    var _f = this.fire;
    this.fire = function () {
      _f.apply(this, arguments);
      if (this.component) this.component.fire.apply(this.component, arguments);
    };
    this.detached = false;
    this.id = params.id;
    this._jsPlumb.div = null;
    this._jsPlumb.initialised = false;
    this._jsPlumb.component = params.component;
    this._jsPlumb.cachedDimensions = null;
    this._jsPlumb.create = params.create;
    this._jsPlumb.initiallyInvisible = params.visible === false;
    this.getElement = function () {
      if (this._jsPlumb.div == null) {
        var div = (this._jsPlumb.div = _jp.getElement(
          this._jsPlumb.create(this._jsPlumb.component)
        ));
        div.style.position = "absolute";
        jsPlumb.addClass(
          div,
          this._jsPlumb.instance.overlayClass +
            " " +
            (this.cssClass
              ? this.cssClass
              : params.cssClass
              ? params.cssClass
              : "")
        );
        this._jsPlumb.instance.appendElement(div);
        this._jsPlumb.instance.getId(div);
        this.canvas = div;
        var ts = "translate(-50%, -50%)";
        div.style.webkitTransform = ts;
        div.style.mozTransform = ts;
        div.style.msTransform = ts;
        div.style.oTransform = ts;
        div.style.transform = ts;
        div._jsPlumb = this;
        if (params.visible === false) div.style.display = "none";
      }
      return this._jsPlumb.div;
    };
    this.draw = function (
      component,
      currentConnectionPaintStyle,
      absolutePosition
    ) {
      var td = _getDimensions(this);
      if (td != null && td.length === 2) {
        var cxy = { x: 0, y: 0 };
        if (absolutePosition)
          cxy = { x: absolutePosition[0], y: absolutePosition[1] };
        else if (component.pointOnPath) {
          var loc = this.loc;
          var absolute = false;
          if (_ju.isString(this.loc) || this.loc < 0 || this.loc > 1) {
            loc = parseInt(this.loc, 10);
            absolute = true;
          }
          cxy = component.pointOnPath(loc, absolute);
        } else {
          var locToUse =
            this.loc.constructor === Array ? this.loc : this.endpointLoc;
          cxy = { x: locToUse[0] * component.w, y: locToUse[1] * component.h };
        }
        var minx = cxy.x - td[0] / 2;
        var miny = cxy.y - td[1] / 2;
        return {
          component,
          d: { minx, miny, td, cxy },
          minX: minx,
          maxX: minx + td[0],
          minY: miny,
          maxY: miny + td[1],
        };
      } else return { minX: 0, maxX: 0, minY: 0, maxY: 0 };
    };
  };
  _ju.extend(AbstractDOMOverlay, [_jp.jsPlumbUIComponent, AbstractOverlay], {
    getDimensions: function () {
      return [1, 1];
    },
    setVisible: function (state) {
      if (this._jsPlumb.div) {
        this._jsPlumb.div.style.display = state ? "block" : "none";
        if (state && this._jsPlumb.initiallyInvisible) {
          _getDimensions(this, true);
          this.component.repaint();
          this._jsPlumb.initiallyInvisible = false;
        }
      }
    },
    clearCachedDimensions: function () {
      this._jsPlumb.cachedDimensions = null;
    },
    cleanup: function (force) {
      if (force) {
        if (this._jsPlumb.div != null) {
          this._jsPlumb.div._jsPlumb = null;
          this._jsPlumb.instance.removeElement(this._jsPlumb.div);
        }
      } else {
        if (this._jsPlumb && this._jsPlumb.div && this._jsPlumb.div.parentNode)
          this._jsPlumb.div.parentNode.removeChild(this._jsPlumb.div);
        this.detached = true;
      }
    },
    reattach: function (instance, component) {
      if (this._jsPlumb.div != null)
        instance.getContainer().appendChild(this._jsPlumb.div);
      this.detached = false;
    },
    computeMaxSize: function () {
      var td = _getDimensions(this);
      return Math.max(td[0], td[1]);
    },
    paint: function (p, containerExtents) {
      if (!this._jsPlumb.initialised) {
        this.getElement();
        p.component.appendDisplayElement(this._jsPlumb.div);
        this._jsPlumb.initialised = true;
        if (this.detached)
          this._jsPlumb.div.parentNode.removeChild(this._jsPlumb.div);
      }
      this._jsPlumb.div.style.left = p.component.x + p.d.minx + "px";
      this._jsPlumb.div.style.top = p.component.y + p.d.miny + "px";
    },
  });
  _jp.Overlays.Custom = function (params) {
    this.type = "Custom";
    AbstractDOMOverlay.apply(this, arguments);
  };
  _ju.extend(_jp.Overlays.Custom, AbstractDOMOverlay);
  _jp.Overlays.GuideLines = function () {
    var self = this;
    self.length = 50;
    self.strokeWidth = 5;
    this.type = "GuideLines";
    AbstractOverlay.apply(this, arguments);
    _jp.jsPlumbUIComponent.apply(this, arguments);
    this.draw = function (connector, currentConnectionPaintStyle) {
      var head = connector.pointAlongPathFrom(self.loc, self.length / 2);
      var mid = connector.pointOnPath(self.loc);
      var tail = _jg.pointOnLine(head, mid, self.length);
      var tailLine = _jg.perpendicularLineTo(head, tail, 40);
      var headLine = _jg.perpendicularLineTo(tail, head, 20);
      return {
        connector,
        head,
        tail,
        headLine,
        tailLine,
        minX: Math.min(head.x, tail.x, headLine[0].x, headLine[1].x),
        minY: Math.min(head.y, tail.y, headLine[0].y, headLine[1].y),
        maxX: Math.max(head.x, tail.x, headLine[0].x, headLine[1].x),
        maxY: Math.max(head.y, tail.y, headLine[0].y, headLine[1].y),
      };
    };
  };
  _jp.Overlays.Label = function (params) {
    this.labelStyle = params.labelStyle;
    var labelWidth = null;
    var labelHeight = null;
    var labelText = null;
    var labelPadding = null;
    this.cssClass = this.labelStyle != null ? this.labelStyle.cssClass : null;
    var p = _jp.extend(
      {
        create: function () {
          return _jp.createElement("div");
        },
      },
      params
    );
    _jp.Overlays.Custom.call(this, p);
    this.type = "Label";
    this.label = params.label || "";
    this.labelText = null;
    if (this.labelStyle) {
      var el = this.getElement();
      this.labelStyle.font = this.labelStyle.font || "12px sans-serif";
      el.style.font = this.labelStyle.font;
      el.style.color = this.labelStyle.color || "black";
      if (this.labelStyle.fill) el.style.background = this.labelStyle.fill;
      if (this.labelStyle.borderWidth > 0) {
        var dStyle = this.labelStyle.borderStyle
          ? this.labelStyle.borderStyle
          : "black";
        el.style.border = this.labelStyle.borderWidth + "px solid " + dStyle;
      }
      if (this.labelStyle.padding) el.style.padding = this.labelStyle.padding;
    }
  };
  _ju.extend(_jp.Overlays.Label, _jp.Overlays.Custom, {
    cleanup: function (force) {
      if (force) {
        this.div = null;
        this.label = null;
        this.labelText = null;
        this.cssClass = null;
        this.labelStyle = null;
      }
    },
    getLabel: function () {
      return this.label;
    },
    setLabel: function (l) {
      this.label = l;
      this.labelText = null;
      this.clearCachedDimensions();
      this.update();
      this.component.repaint();
    },
    getDimensions: function () {
      this.update();
      return AbstractDOMOverlay.prototype.getDimensions.apply(this, arguments);
    },
    update: function () {
      if (typeof this.label === "function") {
        var lt = this.label(this);
        this.getElement().innerHTML = lt.replace(/\r\n/g, "\x3cbr/\x3e");
      } else if (this.labelText == null) {
        this.labelText = this.label;
        this.getElement().innerHTML = this.labelText.replace(
          /\r\n/g,
          "\x3cbr/\x3e"
        );
      }
    },
    updateFrom: function (d) {
      if (d.label != null) this.setLabel(d.label);
    },
  });
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _ju = root.jsPlumbUtil;
  var _jpi = root.jsPlumbInstance;
  var GROUP_COLLAPSED_CLASS = "jtk-group-collapsed";
  var GROUP_EXPANDED_CLASS = "jtk-group-expanded";
  var GROUP_CONTAINER_SELECTOR = "[jtk-group-content]";
  var ELEMENT_DRAGGABLE_EVENT = "elementDraggable";
  var STOP = "stop";
  var REVERT = "revert";
  var GROUP_MANAGER = "_groupManager";
  var GROUP = "_jsPlumbGroup";
  var GROUP_DRAG_SCOPE = "_jsPlumbGroupDrag";
  var EVT_CHILD_ADDED = "group:addMember";
  var EVT_CHILD_REMOVED = "group:removeMember";
  var EVT_GROUP_ADDED = "group:add";
  var EVT_GROUP_REMOVED = "group:remove";
  var EVT_EXPAND = "group:expand";
  var EVT_COLLAPSE = "group:collapse";
  var EVT_GROUP_DRAG_STOP = "groupDragStop";
  var EVT_CONNECTION_MOVED = "connectionMoved";
  var EVT_INTERNAL_CONNECTION_DETACHED = "internal.connectionDetached";
  var CMD_REMOVE_ALL = "removeAll";
  var CMD_ORPHAN_ALL = "orphanAll";
  var CMD_SHOW = "show";
  var CMD_HIDE = "hide";
  var GroupManager = function (_jsPlumb) {
    function isDescendant(el, parentEl) {
      var c = _jsPlumb.getContainer();
      var abort = false;
      var g = null;
      for (var child = null; !abort; )
        if (el == null || el === c) return false;
        else if (el === parentEl) return true;
        else el = el.parentNode;
    }
    function _cleanupDetachedConnection(conn) {
      delete conn.proxies;
      var group = _connectionSourceMap[conn.id];
      if (group != null) {
        var f = function (c) {
          return c.id === conn.id;
        };
        _ju.removeWithFunction(group.connections.source, f);
        _ju.removeWithFunction(group.connections.target, f);
        delete _connectionSourceMap[conn.id];
      }
      group = _connectionTargetMap[conn.id];
      if (group != null) {
        f = function (c) {
          return c.id === conn.id;
        };
        _ju.removeWithFunction(group.connections.source, f);
        _ju.removeWithFunction(group.connections.target, f);
        delete _connectionTargetMap[conn.id];
      }
    }
    function _setVisible(group, state) {
      var m = group.getEl().querySelectorAll(".jtk-managed");
      for (var i = 0; i < m.length; i++)
        _jsPlumb[state ? CMD_SHOW : CMD_HIDE](m[i], true);
    }
    function _updateConnectionsForGroup(group) {
      var members = group.getMembers().slice();
      var childMembers = [];
      for (var i$jscomp$0 = 0; i$jscomp$0 < members.length; i$jscomp$0++)
        Array.prototype.push.apply(
          childMembers,
          members[i$jscomp$0].querySelectorAll(".jtk-managed")
        );
      Array.prototype.push.apply(members, childMembers);
      var c1 = _jsPlumb.getConnections({ source: members, scope: "*" }, true);
      var c2 = _jsPlumb.getConnections({ target: members, scope: "*" }, true);
      var processed = {};
      group.connections.source.length = 0;
      group.connections.target.length = 0;
      var oneSet = function (c) {
        for (var i = 0; i < c.length; i++) {
          if (processed[c[i].id]) continue;
          processed[c[i].id] = true;
          var gs = _jsPlumb.getGroupFor(c[i].source);
          var gt = _jsPlumb.getGroupFor(c[i].target);
          if (gs === group) {
            if (gt !== group) group.connections.source.push(c[i]);
            _connectionSourceMap[c[i].id] = group;
          } else if (gt === group) {
            group.connections.target.push(c[i]);
            _connectionTargetMap[c[i].id] = group;
          }
        }
      };
      oneSet(c1);
      oneSet(c2);
    }
    var _managedGroups = {};
    var _connectionSourceMap = {};
    var _connectionTargetMap = {};
    var self = this;
    _jsPlumb.bind("connection", function (p) {
      var sourceGroup = _jsPlumb.getGroupFor(p.source);
      var targetGroup = _jsPlumb.getGroupFor(p.target);
      if (
        sourceGroup != null &&
        targetGroup != null &&
        sourceGroup === targetGroup
      ) {
        _connectionSourceMap[p.connection.id] = sourceGroup;
        _connectionTargetMap[p.connection.id] = sourceGroup;
      } else {
        if (sourceGroup != null) {
          _ju.suggest(sourceGroup.connections.source, p.connection);
          _connectionSourceMap[p.connection.id] = sourceGroup;
        }
        if (targetGroup != null) {
          _ju.suggest(targetGroup.connections.target, p.connection);
          _connectionTargetMap[p.connection.id] = targetGroup;
        }
      }
    });
    _jsPlumb.bind(EVT_INTERNAL_CONNECTION_DETACHED, function (p) {
      _cleanupDetachedConnection(p.connection);
    });
    _jsPlumb.bind(EVT_CONNECTION_MOVED, function (p) {
      var connMap = p.index === 0 ? _connectionSourceMap : _connectionTargetMap;
      var group = connMap[p.connection.id];
      if (group) {
        var list = group.connections[p.index === 0 ? "source" : "target"];
        var idx = list.indexOf(p.connection);
        if (idx !== -1) list.splice(idx, 1);
      }
    });
    this.addGroup = function (group) {
      _jsPlumb.addClass(group.getEl(), GROUP_EXPANDED_CLASS);
      _managedGroups[group.id] = group;
      group.manager = this;
      _updateConnectionsForGroup(group);
      _jsPlumb.fire(EVT_GROUP_ADDED, { group });
    };
    this.addToGroup = function (group, el, doNotFireEvent) {
      group = this.getGroup(group);
      if (group) {
        var groupEl = group.getEl();
        if (el._isJsPlumbGroup) return;
        var currentGroup = el._jsPlumbGroup;
        if (currentGroup !== group) {
          _jsPlumb.removeFromDragSelection(el);
          var elpos = _jsPlumb.getOffset(el, true);
          var cpos = group.collapsed
            ? _jsPlumb.getOffset(groupEl, true)
            : _jsPlumb.getOffset(group.getDragArea(), true);
          if (currentGroup != null) {
            currentGroup.remove(el, false, doNotFireEvent, false, group);
            self.updateConnectionsForGroup(currentGroup);
          }
          group.add(el, doNotFireEvent);
          var handleDroppedConnections = function (list, index) {
            var oidx = index === 0 ? 1 : 0;
            list.each(function (c) {
              c.setVisible(false);
              if (c.endpoints[oidx].element._jsPlumbGroup === group) {
                c.endpoints[oidx].setVisible(false);
                _expandConnection(c, oidx, group);
              } else {
                c.endpoints[index].setVisible(false);
                _collapseConnection(c, index, group);
              }
            });
          };
          if (group.collapsed) {
            handleDroppedConnections(_jsPlumb.select({ source: el }), 0);
            handleDroppedConnections(_jsPlumb.select({ target: el }), 1);
          }
          var elId = _jsPlumb.getId(el);
          _jsPlumb.dragManager.setParent(
            el,
            elId,
            groupEl,
            _jsPlumb.getId(groupEl),
            elpos
          );
          var newPosition = {
            left: elpos.left - cpos.left,
            top: elpos.top - cpos.top,
          };
          _jsPlumb.setPosition(el, newPosition);
          _jsPlumb.dragManager.revalidateParent(el, elId, elpos);
          self.updateConnectionsForGroup(group);
          _jsPlumb.revalidate(elId);
          if (!doNotFireEvent) {
            var p = { group, el, pos: newPosition };
            if (currentGroup) p.sourceGroup = currentGroup;
            _jsPlumb.fire(EVT_CHILD_ADDED, p);
          }
        }
      }
    };
    this.removeFromGroup = function (group, el, doNotFireEvent) {
      group = this.getGroup(group);
      if (group) {
        if (group.collapsed) {
          var _expandSet = function (conns, index) {
            for (var i = 0; i < conns.length; i++) {
              var c = conns[i];
              if (c.proxies)
                for (var j = 0; j < c.proxies.length; j++)
                  if (c.proxies[j] != null) {
                    var proxiedElement = c.proxies[j].originalEp.element;
                    if (
                      proxiedElement === el ||
                      isDescendant(proxiedElement, el)
                    )
                      _expandConnection(c, index, group);
                  }
            }
          };
          _expandSet(group.connections.source.slice(), 0);
          _expandSet(group.connections.target.slice(), 1);
        }
        group.remove(el, null, doNotFireEvent);
      }
    };
    this.getGroup = function (groupId) {
      var group = groupId;
      if (_ju.isString(groupId)) {
        group = _managedGroups[groupId];
        if (group == null)
          throw new TypeError("No such group [" + groupId + "]");
      }
      return group;
    };
    this.getGroups = function () {
      var o = [];
      for (var g in _managedGroups) o.push(_managedGroups[g]);
      return o;
    };
    this.removeGroup = function (
      group,
      deleteMembers,
      manipulateDOM,
      doNotFireEvent
    ) {
      group = this.getGroup(group);
      this.expandGroup(group, true);
      var newPositions = group[deleteMembers ? CMD_REMOVE_ALL : CMD_ORPHAN_ALL](
        manipulateDOM,
        doNotFireEvent
      );
      _jsPlumb.remove(group.getEl());
      delete _managedGroups[group.id];
      delete _jsPlumb._groups[group.id];
      _jsPlumb.fire(EVT_GROUP_REMOVED, { group });
      return newPositions;
    };
    this.removeAllGroups = function (
      deleteMembers,
      manipulateDOM,
      doNotFireEvent
    ) {
      for (var g in _managedGroups)
        this.removeGroup(
          _managedGroups[g],
          deleteMembers,
          manipulateDOM,
          doNotFireEvent
        );
    };
    var _collapseConnection = function (c$jscomp$0, index$jscomp$0, group) {
      var otherEl = c$jscomp$0.endpoints[index$jscomp$0 === 0 ? 1 : 0].element;
      if (
        otherEl[GROUP] &&
        !otherEl[GROUP].shouldProxy() &&
        otherEl[GROUP].collapsed
      )
        return;
      var groupEl = group.getEl();
      var groupElId = _jsPlumb.getId(groupEl);
      _jsPlumb.proxyConnection(
        c$jscomp$0,
        index$jscomp$0,
        groupEl,
        groupElId,
        function (c, index) {
          return group.getEndpoint(c, index);
        },
        function (c, index) {
          return group.getAnchor(c, index);
        }
      );
    };
    this.collapseGroup = function (group) {
      group = this.getGroup(group);
      if (group == null || group.collapsed) return;
      var groupEl = group.getEl();
      _setVisible(group, false);
      if (group.shouldProxy()) {
        var _collapseSet = function (conns, index) {
          for (var i = 0; i < conns.length; i++) {
            var c = conns[i];
            _collapseConnection(c, index, group);
          }
        };
        _collapseSet(group.connections.source, 0);
        _collapseSet(group.connections.target, 1);
      }
      group.collapsed = true;
      _jsPlumb.removeClass(groupEl, GROUP_EXPANDED_CLASS);
      _jsPlumb.addClass(groupEl, GROUP_COLLAPSED_CLASS);
      _jsPlumb.revalidate(groupEl);
      _jsPlumb.fire(EVT_COLLAPSE, { group });
    };
    var _expandConnection = function (c, index, group) {
      _jsPlumb.unproxyConnection(c, index, _jsPlumb.getId(group.getEl()));
    };
    this.expandGroup = function (group, doNotFireEvent) {
      group = this.getGroup(group);
      if (group == null || !group.collapsed) return;
      var groupEl = group.getEl();
      _setVisible(group, true);
      if (group.shouldProxy()) {
        var _expandSet = function (conns, index) {
          for (var i = 0; i < conns.length; i++) {
            var c = conns[i];
            _expandConnection(c, index, group);
          }
        };
        _expandSet(group.connections.source, 0);
        _expandSet(group.connections.target, 1);
      }
      group.collapsed = false;
      _jsPlumb.addClass(groupEl, GROUP_EXPANDED_CLASS);
      _jsPlumb.removeClass(groupEl, GROUP_COLLAPSED_CLASS);
      _jsPlumb.revalidate(groupEl);
      this.repaintGroup(group);
      if (!doNotFireEvent) _jsPlumb.fire(EVT_EXPAND, { group });
    };
    this.repaintGroup = function (group) {
      group = this.getGroup(group);
      var m = group.getMembers();
      for (var i = 0; i < m.length; i++) _jsPlumb.revalidate(m[i]);
    };
    this.updateConnectionsForGroup = _updateConnectionsForGroup;
    this.refreshAllGroups = function () {
      for (var g in _managedGroups) {
        _updateConnectionsForGroup(_managedGroups[g]);
        _jsPlumb.dragManager.updateOffsets(
          _jsPlumb.getId(_managedGroups[g].getEl())
        );
      }
    };
  };
  var Group = function (_jsPlumb, params$jscomp$0) {
    function _findParent(_el) {
      return _el.offsetParent;
    }
    function _isInsideParent(_el, pos) {
      var p = _findParent(_el);
      var s = _jsPlumb.getSize(p);
      var ss = _jsPlumb.getSize(_el);
      var leftEdge = pos[0];
      var rightEdge = leftEdge + ss[0];
      var topEdge = pos[1];
      var bottomEdge = topEdge + ss[1];
      return (
        rightEdge > 0 && leftEdge < s[0] && bottomEdge > 0 && topEdge < s[1]
      );
    }
    function _orphan(_el) {
      var id = _jsPlumb.getId(_el);
      var pos = _jsPlumb.getOffset(_el);
      _el.parentNode.removeChild(_el);
      _jsPlumb.getContainer().appendChild(_el);
      _jsPlumb.setPosition(_el, pos);
      _unbindDragHandlers(_el);
      _jsPlumb.dragManager.clearParent(_el, id);
      return [id, pos];
    }
    function _pruneOrOrphan(p) {
      function _one(el, left, top) {
        var orphanedPosition = null;
        if (!_isInsideParent(el, [left, top])) {
          var group = el._jsPlumbGroup;
          if (prune) _jsPlumb.remove(el);
          else orphanedPosition = _orphan(el);
          group.remove(el);
        }
        return orphanedPosition;
      }
      var out = [];
      for (var i = 0; i < p.selection.length; i++)
        out.push(
          _one(p.selection[i][0], p.selection[i][1].left, p.selection[i][1].top)
        );
      return out.length === 1 ? out[0] : out;
    }
    function _revalidate(_el) {
      var id = _jsPlumb.getId(_el);
      _jsPlumb.revalidate(_el);
      _jsPlumb.dragManager.revalidateParent(_el, id);
    }
    function _unbindDragHandlers(_el) {
      if (!_el._katavorioDrag) return;
      if (prune || orphan) _el._katavorioDrag.off(STOP, _pruneOrOrphan);
      if (!prune && !orphan && revert) {
        _el._katavorioDrag.off(REVERT, _revalidate);
        _el._katavorioDrag.setRevert(null);
      }
    }
    function _bindDragHandlers(_el) {
      if (!_el._katavorioDrag) return;
      if (prune || orphan) _el._katavorioDrag.on(STOP, _pruneOrOrphan);
      if (constrain) _el._katavorioDrag.setConstrain(true);
      if (ghost) _el._katavorioDrag.setUseGhostProxy(true);
      if (!prune && !orphan && revert) {
        _el._katavorioDrag.on(REVERT, _revalidate);
        _el._katavorioDrag.setRevert(function (__el, pos) {
          return !_isInsideParent(__el, pos);
        });
      }
    }
    var self = this;
    var el$jscomp$0 = params$jscomp$0.el;
    this.getEl = function () {
      return el$jscomp$0;
    };
    this.id = params$jscomp$0.id || _ju.uuid();
    el$jscomp$0._isJsPlumbGroup = true;
    var getDragArea = (this.getDragArea = function () {
      var da = _jsPlumb.getSelector(el$jscomp$0, GROUP_CONTAINER_SELECTOR);
      return da && da.length > 0 ? da[0] : el$jscomp$0;
    });
    var ghost = params$jscomp$0.ghost === true;
    var constrain = ghost || params$jscomp$0.constrain === true;
    var revert = params$jscomp$0.revert !== false;
    var orphan = params$jscomp$0.orphan === true;
    var prune = params$jscomp$0.prune === true;
    var dropOverride = params$jscomp$0.dropOverride === true;
    var proxied = params$jscomp$0.proxied !== false;
    var elements = [];
    this.connections = { source: [], target: [], internal: [] };
    this.getAnchor = function (conn, endpointIndex) {
      return params$jscomp$0.anchor || "Continuous";
    };
    this.getEndpoint = function (conn, endpointIndex) {
      return params$jscomp$0.endpoint || ["Dot", { radius: 10 }];
    };
    this.collapsed = false;
    if (params$jscomp$0.draggable !== false) {
      var opts = {
        drag: function () {
          for (var i = 0; i < elements.length; i++) _jsPlumb.draw(elements[i]);
        },
        stop: function (params) {
          _jsPlumb.fire(
            EVT_GROUP_DRAG_STOP,
            jsPlumb.extend(params, { group: self })
          );
        },
        scope: GROUP_DRAG_SCOPE,
      };
      if (params$jscomp$0.dragOptions)
        root.jsPlumb.extend(opts, params$jscomp$0.dragOptions);
      _jsPlumb.draggable(params$jscomp$0.el, opts);
    }
    if (params$jscomp$0.droppable !== false)
      _jsPlumb.droppable(params$jscomp$0.el, {
        drop: function (p) {
          var el = p.drag.el;
          if (el._isJsPlumbGroup) return;
          var currentGroup = el._jsPlumbGroup;
          if (currentGroup !== self) {
            if (currentGroup != null)
              if (currentGroup.overrideDrop(el, self)) return;
            _jsPlumb.getGroupManager().addToGroup(self, el, false);
          }
        },
      });
    var _each = function (_el, fn) {
      var els = _el.nodeType == null ? _el : [_el];
      for (var i = 0; i < els.length; i++) fn(els[i]);
    };
    this.overrideDrop = function (_el, targetGroup) {
      return dropOverride && (revert || prune || orphan);
    };
    this.add = function (_el, doNotFireEvent) {
      var dragArea = getDragArea();
      _each(_el, function (__el) {
        if (__el._jsPlumbGroup != null)
          if (__el._jsPlumbGroup === self) return;
          else __el._jsPlumbGroup.remove(__el, true, doNotFireEvent, false);
        __el._jsPlumbGroup = self;
        elements.push(__el);
        if (_jsPlumb.isAlreadyDraggable(__el)) _bindDragHandlers(__el);
        if (__el.parentNode !== dragArea) dragArea.appendChild(__el);
      });
      _jsPlumb.getGroupManager().updateConnectionsForGroup(self);
    };
    this.remove = function (
      el,
      manipulateDOM,
      doNotFireEvent,
      doNotUpdateConnections,
      targetGroup
    ) {
      _each(el, function (__el) {
        if (__el._jsPlumbGroup === self) {
          delete __el._jsPlumbGroup;
          _ju.removeWithFunction(elements, function (e) {
            return e === __el;
          });
          if (manipulateDOM)
            try {
              self.getDragArea().removeChild(__el);
            } catch (e) {
              jsPlumbUtil.log("Could not remove element from Group " + e);
            }
          _unbindDragHandlers(__el);
          if (!doNotFireEvent) {
            var p = { group: self, el: __el };
            if (targetGroup) p.targetGroup = targetGroup;
            _jsPlumb.fire(EVT_CHILD_REMOVED, p);
          }
        }
      });
      if (!doNotUpdateConnections)
        _jsPlumb.getGroupManager().updateConnectionsForGroup(self);
    };
    this.removeAll = function (manipulateDOM, doNotFireEvent) {
      var i = 0;
      for (var l = elements.length; i < l; i++) {
        var el = elements[0];
        self.remove(el, manipulateDOM, doNotFireEvent, true);
        _jsPlumb.remove(el, true);
      }
      elements.length = 0;
      _jsPlumb.getGroupManager().updateConnectionsForGroup(self);
    };
    this.orphanAll = function () {
      var orphanedPositions = {};
      for (var i = 0; i < elements.length; i++) {
        var newPosition = _orphan(elements[i]);
        orphanedPositions[newPosition[0]] = newPosition[1];
      }
      elements.length = 0;
      return orphanedPositions;
    };
    this.getMembers = function () {
      return elements;
    };
    el$jscomp$0[GROUP] = this;
    _jsPlumb.bind(
      ELEMENT_DRAGGABLE_EVENT,
      function (dragParams) {
        if (dragParams.el._jsPlumbGroup === this)
          _bindDragHandlers(dragParams.el);
      }.bind(this)
    );
    this.shouldProxy = function () {
      return proxied;
    };
    _jsPlumb.getGroupManager().addGroup(this);
  };
  _jpi.prototype.addGroup = function (params) {
    var j = this;
    j._groups = j._groups || {};
    if (j._groups[params.id] != null)
      throw new TypeError(
        "cannot create Group [" + params.id + "]; a Group with that ID exists"
      );
    if (params.el[GROUP] != null)
      throw new TypeError(
        "cannot create Group [" +
          params.id +
          "]; the given element is already a Group"
      );
    var group = new Group(j, params);
    j._groups[group.id] = group;
    if (params.collapsed) this.collapseGroup(group);
    return group;
  };
  _jpi.prototype.addToGroup = function (group, el, doNotFireEvent) {
    var _one = function (_el) {
      var id = this.getId(_el);
      this.manage(id, _el);
      this.getGroupManager().addToGroup(group, _el, doNotFireEvent);
    }.bind(this);
    if (Array.isArray(el)) for (var i = 0; i < el.length; i++) _one(el[i]);
    else _one(el);
  };
  _jpi.prototype.removeFromGroup = function (group, el, doNotFireEvent) {
    this.getGroupManager().removeFromGroup(group, el, doNotFireEvent);
    this.getContainer().appendChild(el);
  };
  _jpi.prototype.removeGroup = function (
    group,
    deleteMembers,
    manipulateDOM,
    doNotFireEvent
  ) {
    return this.getGroupManager().removeGroup(
      group,
      deleteMembers,
      manipulateDOM,
      doNotFireEvent
    );
  };
  _jpi.prototype.removeAllGroups = function (
    deleteMembers,
    manipulateDOM,
    doNotFireEvent
  ) {
    this.getGroupManager().removeAllGroups(
      deleteMembers,
      manipulateDOM,
      doNotFireEvent
    );
  };
  _jpi.prototype.getGroup = function (groupId) {
    return this.getGroupManager().getGroup(groupId);
  };
  _jpi.prototype.getGroups = function () {
    return this.getGroupManager().getGroups();
  };
  _jpi.prototype.expandGroup = function (group) {
    this.getGroupManager().expandGroup(group);
  };
  _jpi.prototype.collapseGroup = function (groupId) {
    this.getGroupManager().collapseGroup(groupId);
  };
  _jpi.prototype.repaintGroup = function (group) {
    this.getGroupManager().repaintGroup(group);
  };
  _jpi.prototype.toggleGroup = function (group) {
    group = this.getGroupManager().getGroup(group);
    if (group != null)
      this.getGroupManager()[group.collapsed ? "expandGroup" : "collapseGroup"](
        group
      );
  };
  _jpi.prototype.getGroupManager = function () {
    var mgr = this[GROUP_MANAGER];
    if (mgr == null) mgr = this[GROUP_MANAGER] = new GroupManager(this);
    return mgr;
  };
  _jpi.prototype.removeGroupManager = function () {
    delete this[GROUP_MANAGER];
  };
  _jpi.prototype.getGroupFor = function (el) {
    el = this.getElement(el);
    if (el) {
      var c = this.getContainer();
      var abort = false;
      var g = null;
      for (var child = null; !abort; )
        if (el == null || el === c) abort = true;
        else if (el[GROUP]) {
          g = el[GROUP];
          child = el;
          abort = true;
        } else el = el.parentNode;
      return g;
    }
  };
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var STRAIGHT = "Straight";
  var ARC = "Arc";
  var Flowchart = function (params$jscomp$0) {
    this.type = "Flowchart";
    params$jscomp$0 = params$jscomp$0 || {};
    params$jscomp$0.stub =
      params$jscomp$0.stub == null ? 30 : params$jscomp$0.stub;
    var segments$jscomp$0;
    var _super = _jp.Connectors.AbstractConnector.apply(this, arguments);
    var midpoint =
      params$jscomp$0.midpoint == null ? 0.5 : params$jscomp$0.midpoint;
    var alwaysRespectStubs = params$jscomp$0.alwaysRespectStubs === true;
    var lastx = null;
    var lasty = null;
    var lastOrientation;
    var cornerRadius =
      params$jscomp$0.cornerRadius != null ? params$jscomp$0.cornerRadius : 0;
    var loopbackRadius = params$jscomp$0.loopbackRadius || 25;
    var isLoopbackCurrently = false;
    var sgn = function (n) {
      return n < 0 ? -1 : n === 0 ? 0 : 1;
    };
    var segmentDirections = function (segment) {
      return [sgn(segment[2] - segment[0]), sgn(segment[3] - segment[1])];
    };
    var addSegment = function (segments, x, y, paintInfo) {
      if (lastx === x && lasty === y) return;
      var lx = lastx == null ? paintInfo.sx : lastx;
      var ly = lasty == null ? paintInfo.sy : lasty;
      var o = lx === x ? "v" : "h";
      lastx = x;
      lasty = y;
      segments.push([lx, ly, x, y, o]);
    };
    var segLength = function (s) {
      return Math.sqrt(Math.pow(s[0] - s[2], 2) + Math.pow(s[1] - s[3], 2));
    };
    var _cloneArray = function (a) {
      var _a = [];
      _a.push.apply(_a, a);
      return _a;
    };
    var writeSegments = function (conn, segments, paintInfo) {
      var current = null;
      for (var i = 0; i < segments.length - 1; i++) {
        current = current || _cloneArray(segments[i]);
        var next = _cloneArray(segments[i + 1]);
        var currentDirection = segmentDirections(current);
        var nextDirection = segmentDirections(next);
        if (cornerRadius > 0 && current[4] !== next[4]) {
          var minSegLength = Math.min(segLength(current), segLength(next));
          var radiusToUse = Math.min(cornerRadius, minSegLength / 2);
          current[2] -= currentDirection[0] * radiusToUse;
          current[3] -= currentDirection[1] * radiusToUse;
          next[0] += nextDirection[0] * radiusToUse;
          next[1] += nextDirection[1] * radiusToUse;
          var ac =
            (currentDirection[1] === nextDirection[0] &&
              nextDirection[0] === 1) ||
            (currentDirection[1] === nextDirection[0] &&
              nextDirection[0] === 0 &&
              currentDirection[0] !== nextDirection[1]) ||
            (currentDirection[1] === nextDirection[0] &&
              nextDirection[0] === -1);
          var sgny = next[1] > current[3] ? 1 : -1;
          var sgnx = next[0] > current[2] ? 1 : -1;
          var sgnEqual = sgny === sgnx;
          var cx =
            (sgnEqual && ac) || (!sgnEqual && !ac) ? next[0] : current[2];
          var cy =
            (sgnEqual && ac) || (!sgnEqual && !ac) ? current[3] : next[1];
          _super.addSegment(conn, STRAIGHT, {
            x1: current[0],
            y1: current[1],
            x2: current[2],
            y2: current[3],
          });
          _super.addSegment(conn, ARC, {
            r: radiusToUse,
            x1: current[2],
            y1: current[3],
            x2: next[0],
            y2: next[1],
            cx,
            cy,
            ac,
          });
        } else {
          var dx =
            current[2] === current[0]
              ? 0
              : current[2] > current[0]
              ? paintInfo.lw / 2
              : -(paintInfo.lw / 2);
          var dy =
            current[3] === current[1]
              ? 0
              : current[3] > current[1]
              ? paintInfo.lw / 2
              : -(paintInfo.lw / 2);
          _super.addSegment(conn, STRAIGHT, {
            x1: current[0] - dx,
            y1: current[1] - dy,
            x2: current[2] + dx,
            y2: current[3] + dy,
          });
        }
        current = next;
      }
      if (next != null)
        _super.addSegment(conn, STRAIGHT, {
          x1: next[0],
          y1: next[1],
          x2: next[2],
          y2: next[3],
        });
    };
    this._compute = function (paintInfo, params) {
      segments$jscomp$0 = [];
      lastx = null;
      lasty = null;
      lastOrientation = null;
      var commonStubCalculator = function () {
        return [
          paintInfo.startStubX,
          paintInfo.startStubY,
          paintInfo.endStubX,
          paintInfo.endStubY,
        ];
      };
      var stubCalculators = {
        perpendicular: commonStubCalculator,
        orthogonal: commonStubCalculator,
        opposite: function (axis) {
          var pi = paintInfo;
          var idx = axis === "x" ? 0 : 1;
          var areInProximity = {
            x: function () {
              return (
                (pi.so[idx] === 1 &&
                  ((pi.startStubX > pi.endStubX && pi.tx > pi.startStubX) ||
                    (pi.sx > pi.endStubX && pi.tx > pi.sx))) ||
                (pi.so[idx] === -1 &&
                  ((pi.startStubX < pi.endStubX && pi.tx < pi.startStubX) ||
                    (pi.sx < pi.endStubX && pi.tx < pi.sx)))
              );
            },
            y: function () {
              return (
                (pi.so[idx] === 1 &&
                  ((pi.startStubY > pi.endStubY && pi.ty > pi.startStubY) ||
                    (pi.sy > pi.endStubY && pi.ty > pi.sy))) ||
                (pi.so[idx] === -1 &&
                  ((pi.startStubY < pi.endStubY && pi.ty < pi.startStubY) ||
                    (pi.sy < pi.endStubY && pi.ty < pi.sy)))
              );
            },
          };
          if (!alwaysRespectStubs && areInProximity[axis]())
            return {
              x: [
                (paintInfo.sx + paintInfo.tx) / 2,
                paintInfo.startStubY,
                (paintInfo.sx + paintInfo.tx) / 2,
                paintInfo.endStubY,
              ],
              y: [
                paintInfo.startStubX,
                (paintInfo.sy + paintInfo.ty) / 2,
                paintInfo.endStubX,
                (paintInfo.sy + paintInfo.ty) / 2,
              ],
            }[axis];
          else
            return [
              paintInfo.startStubX,
              paintInfo.startStubY,
              paintInfo.endStubX,
              paintInfo.endStubY,
            ];
        },
      };
      var stubs = stubCalculators[paintInfo.anchorOrientation](
        paintInfo.sourceAxis
      );
      var idx$jscomp$0 = paintInfo.sourceAxis === "x" ? 0 : 1;
      var oidx = paintInfo.sourceAxis === "x" ? 1 : 0;
      var ss$jscomp$0 = stubs[idx$jscomp$0];
      var oss = stubs[oidx];
      var es$jscomp$0 = stubs[idx$jscomp$0 + 2];
      var oes = stubs[oidx + 2];
      addSegment(segments$jscomp$0, stubs[0], stubs[1], paintInfo);
      var midx =
        paintInfo.startStubX +
        (paintInfo.endStubX - paintInfo.startStubX) * midpoint;
      var midy =
        paintInfo.startStubY +
        (paintInfo.endStubY - paintInfo.startStubY) * midpoint;
      var orientations = { x: [0, 1], y: [1, 0] };
      var lineCalculators = {
        perpendicular: function (axis) {
          var pi = paintInfo;
          var sis = {
            x: [
              [[1, 2, 3, 4], null, [2, 1, 4, 3]],
              null,
              [[4, 3, 2, 1], null, [3, 4, 1, 2]],
            ],
            y: [
              [[3, 2, 1, 4], null, [2, 3, 4, 1]],
              null,
              [[4, 1, 2, 3], null, [1, 4, 3, 2]],
            ],
          };
          var stubs = {
            x: [
              [pi.startStubX, pi.endStubX],
              null,
              [pi.endStubX, pi.startStubX],
            ],
            y: [
              [pi.startStubY, pi.endStubY],
              null,
              [pi.endStubY, pi.startStubY],
            ],
          };
          var midLines = {
            x: [
              [midx, pi.startStubY],
              [midx, pi.endStubY],
            ],
            y: [
              [pi.startStubX, midy],
              [pi.endStubX, midy],
            ],
          };
          var linesToEnd = {
            x: [[pi.endStubX, pi.startStubY]],
            y: [[pi.startStubX, pi.endStubY]],
          };
          var startToEnd = {
            x: [
              [pi.startStubX, pi.endStubY],
              [pi.endStubX, pi.endStubY],
            ],
            y: [
              [pi.endStubX, pi.startStubY],
              [pi.endStubX, pi.endStubY],
            ],
          };
          var startToMidToEnd = {
            x: [
              [pi.startStubX, midy],
              [pi.endStubX, midy],
              [pi.endStubX, pi.endStubY],
            ],
            y: [
              [midx, pi.startStubY],
              [midx, pi.endStubY],
              [pi.endStubX, pi.endStubY],
            ],
          };
          var otherStubs = {
            x: [pi.startStubY, pi.endStubY],
            y: [pi.startStubX, pi.endStubX],
          };
          var soIdx = orientations[axis][0];
          var toIdx = orientations[axis][1];
          var _so = pi.so[soIdx] + 1;
          var _to = pi.to[toIdx] + 1;
          var otherFlipped =
            (pi.to[toIdx] === -1 &&
              otherStubs[axis][1] < otherStubs[axis][0]) ||
            (pi.to[toIdx] === 1 && otherStubs[axis][1] > otherStubs[axis][0]);
          var stub1 = stubs[axis][_so][0];
          var stub2 = stubs[axis][_so][1];
          var segmentIndexes = sis[axis][_so][_to];
          if (
            pi.segment === segmentIndexes[3] ||
            (pi.segment === segmentIndexes[2] && otherFlipped)
          )
            return midLines[axis];
          else if (pi.segment === segmentIndexes[2] && stub2 < stub1)
            return linesToEnd[axis];
          else if (
            (pi.segment === segmentIndexes[2] && stub2 >= stub1) ||
            (pi.segment === segmentIndexes[1] && !otherFlipped)
          )
            return startToMidToEnd[axis];
          else if (
            pi.segment === segmentIndexes[0] ||
            (pi.segment === segmentIndexes[1] && otherFlipped)
          )
            return startToEnd[axis];
        },
        orthogonal: function (
          axis,
          startStub,
          otherStartStub,
          endStub,
          otherEndStub
        ) {
          var pi = paintInfo;
          var extent = {
            x:
              pi.so[0] === -1
                ? Math.min(startStub, endStub)
                : Math.max(startStub, endStub),
            y:
              pi.so[1] === -1
                ? Math.min(startStub, endStub)
                : Math.max(startStub, endStub),
          }[axis];
          return {
            x: [
              [extent, otherStartStub],
              [extent, otherEndStub],
              [endStub, otherEndStub],
            ],
            y: [
              [otherStartStub, extent],
              [otherEndStub, extent],
              [otherEndStub, endStub],
            ],
          }[axis];
        },
        opposite: function (axis, ss, oss, es) {
          var pi = paintInfo;
          var otherAxis = { x: "y", y: "x" }[axis];
          var dim = { x: "height", y: "width" }[axis];
          var comparator =
            pi["is" + axis.toUpperCase() + "GreaterThanStubTimes2"];
          if (
            params.sourceEndpoint.elementId === params.targetEndpoint.elementId
          ) {
            var _val =
              oss +
              (1 - params.sourceEndpoint.anchor[otherAxis]) *
                params.sourceInfo[dim] +
              _super.maxStub;
            return {
              x: [
                [ss, _val],
                [es, _val],
              ],
              y: [
                [_val, ss],
                [_val, es],
              ],
            }[axis];
          } else if (
            !comparator ||
            (pi.so[idx$jscomp$0] === 1 && ss > es) ||
            (pi.so[idx$jscomp$0] === -1 && ss < es)
          )
            return {
              x: [
                [ss, midy],
                [es, midy],
              ],
              y: [
                [midx, ss],
                [midx, es],
              ],
            }[axis];
          else if (
            (pi.so[idx$jscomp$0] === 1 && ss < es) ||
            (pi.so[idx$jscomp$0] === -1 && ss > es)
          )
            return {
              x: [
                [midx, pi.sy],
                [midx, pi.ty],
              ],
              y: [
                [pi.sx, midy],
                [pi.tx, midy],
              ],
            }[axis];
        },
      };
      var p = lineCalculators[paintInfo.anchorOrientation](
        paintInfo.sourceAxis,
        ss$jscomp$0,
        oss,
        es$jscomp$0,
        oes
      );
      if (p)
        for (var i = 0; i < p.length; i++)
          addSegment(segments$jscomp$0, p[i][0], p[i][1], paintInfo);
      addSegment(segments$jscomp$0, stubs[2], stubs[3], paintInfo);
      addSegment(segments$jscomp$0, paintInfo.tx, paintInfo.ty, paintInfo);
      writeSegments(this, segments$jscomp$0, paintInfo);
    };
  };
  _jp.Connectors.Flowchart = Flowchart;
  _ju.extend(_jp.Connectors.Flowchart, _jp.Connectors.AbstractConnector);
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  _jp.Connectors.AbstractBezierConnector = function (params) {
    params = params || {};
    var showLoopback = params.showLoopback !== false;
    var curviness = params.curviness || 10;
    var margin = params.margin || 5;
    var proximityLimit = params.proximityLimit || 80;
    var clockwise = params.orientation && params.orientation === "clockwise";
    var loopbackRadius = params.loopbackRadius || 25;
    var isLoopbackCurrently = false;
    this._compute = function (paintInfo, p) {
      var sp = p.sourcePos;
      var tp = p.targetPos;
      var _w = Math.abs(sp[0] - tp[0]);
      var _h = Math.abs(sp[1] - tp[1]);
      if (
        !showLoopback ||
        p.sourceEndpoint.elementId !== p.targetEndpoint.elementId
      ) {
        isLoopbackCurrently = false;
        this._computeBezier(paintInfo, p, sp, tp, _w, _h);
      } else {
        isLoopbackCurrently = true;
        var x1 = p.sourcePos[0];
        var y1 = p.sourcePos[1] - margin;
        var cx = x1;
        var cy = y1 - loopbackRadius;
        var _x = cx - loopbackRadius;
        var _y = cy - loopbackRadius;
        _w = 2 * loopbackRadius;
        _h = 2 * loopbackRadius;
        paintInfo.points[0] = _x;
        paintInfo.points[1] = _y;
        paintInfo.points[2] = _w;
        paintInfo.points[3] = _h;
        _super.addSegment(this, "Arc", {
          loopback: true,
          x1: x1 - _x + 4,
          y1: y1 - _y,
          startAngle: 0,
          endAngle: 2 * Math.PI,
          r: loopbackRadius,
          ac: !clockwise,
          x2: x1 - _x - 4,
          y2: y1 - _y,
          cx: cx - _x,
          cy: cy - _y,
        });
      }
    };
    var _super = _jp.Connectors.AbstractConnector.apply(this, arguments);
    return _super;
  };
  _ju.extend(
    _jp.Connectors.AbstractBezierConnector,
    _jp.Connectors.AbstractConnector
  );
  var Bezier = function (params) {
    params = params || {};
    this.type = "Bezier";
    var _super = _jp.Connectors.AbstractBezierConnector.apply(this, arguments);
    var majorAnchor = params.curviness || 150;
    var minorAnchor = 10;
    this.getCurviness = function () {
      return majorAnchor;
    };
    this._findControlPoint = function (
      point,
      sourceAnchorPosition,
      targetAnchorPosition,
      sourceEndpoint,
      targetEndpoint,
      soo,
      too
    ) {
      var perpendicular = soo[0] !== too[0] || soo[1] === too[1];
      var p = [];
      if (!perpendicular) {
        if (soo[0] === 0)
          p.push(
            sourceAnchorPosition[0] < targetAnchorPosition[0]
              ? point[0] + minorAnchor
              : point[0] - minorAnchor
          );
        else p.push(point[0] - majorAnchor * soo[0]);
        if (soo[1] === 0)
          p.push(
            sourceAnchorPosition[1] < targetAnchorPosition[1]
              ? point[1] + minorAnchor
              : point[1] - minorAnchor
          );
        else p.push(point[1] + majorAnchor * too[1]);
      } else {
        if (too[0] === 0)
          p.push(
            targetAnchorPosition[0] < sourceAnchorPosition[0]
              ? point[0] + minorAnchor
              : point[0] - minorAnchor
          );
        else p.push(point[0] + majorAnchor * too[0]);
        if (too[1] === 0)
          p.push(
            targetAnchorPosition[1] < sourceAnchorPosition[1]
              ? point[1] + minorAnchor
              : point[1] - minorAnchor
          );
        else p.push(point[1] + majorAnchor * soo[1]);
      }
      return p;
    };
    this._computeBezier = function (paintInfo, p, sp, tp, _w, _h) {
      var _sx = sp[0] < tp[0] ? _w : 0;
      var _sy = sp[1] < tp[1] ? _h : 0;
      var _tx = sp[0] < tp[0] ? 0 : _w;
      var _ty = sp[1] < tp[1] ? 0 : _h;
      var _CP = this._findControlPoint(
        [_sx, _sy],
        sp,
        tp,
        p.sourceEndpoint,
        p.targetEndpoint,
        paintInfo.so,
        paintInfo.to
      );
      var _CP2 = this._findControlPoint(
        [_tx, _ty],
        tp,
        sp,
        p.targetEndpoint,
        p.sourceEndpoint,
        paintInfo.to,
        paintInfo.so
      );
      _super.addSegment(this, "Bezier", {
        x1: _sx,
        y1: _sy,
        x2: _tx,
        y2: _ty,
        cp1x: _CP[0],
        cp1y: _CP[1],
        cp2x: _CP2[0],
        cp2y: _CP2[1],
      });
    };
  };
  _jp.Connectors.Bezier = Bezier;
  _ju.extend(Bezier, _jp.Connectors.AbstractBezierConnector);
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var _segment = function (x1, y1, x2, y2) {
    if (x1 <= x2 && y2 <= y1) return 1;
    else if (x1 <= x2 && y1 <= y2) return 2;
    else if (x2 <= x1 && y2 >= y1) return 3;
    return 4;
  };
  var _findControlPoint = function (
    midx,
    midy,
    segment,
    sourceEdge,
    targetEdge,
    dx,
    dy,
    distance,
    proximityLimit
  ) {
    if (distance <= proximityLimit) return [midx, midy];
    if (segment === 1)
      if (sourceEdge[3] <= 0 && targetEdge[3] >= 1)
        return [midx + (sourceEdge[2] < 0.5 ? -1 * dx : dx), midy];
      else if (sourceEdge[2] >= 1 && targetEdge[2] <= 0)
        return [midx, midy + (sourceEdge[3] < 0.5 ? -1 * dy : dy)];
      else return [midx + -1 * dx, midy + -1 * dy];
    else if (segment === 2)
      if (sourceEdge[3] >= 1 && targetEdge[3] <= 0)
        return [midx + (sourceEdge[2] < 0.5 ? -1 * dx : dx), midy];
      else if (sourceEdge[2] >= 1 && targetEdge[2] <= 0)
        return [midx, midy + (sourceEdge[3] < 0.5 ? -1 * dy : dy)];
      else return [midx + dx, midy + -1 * dy];
    else if (segment === 3)
      if (sourceEdge[3] >= 1 && targetEdge[3] <= 0)
        return [midx + (sourceEdge[2] < 0.5 ? -1 * dx : dx), midy];
      else if (sourceEdge[2] <= 0 && targetEdge[2] >= 1)
        return [midx, midy + (sourceEdge[3] < 0.5 ? -1 * dy : dy)];
      else return [midx + -1 * dx, midy + -1 * dy];
    else if (segment === 4)
      if (sourceEdge[3] <= 0 && targetEdge[3] >= 1)
        return [midx + (sourceEdge[2] < 0.5 ? -1 * dx : dx), midy];
      else if (sourceEdge[2] <= 0 && targetEdge[2] >= 1)
        return [midx, midy + (sourceEdge[3] < 0.5 ? -1 * dy : dy)];
      else return [midx + dx, midy + -1 * dy];
  };
  var StateMachine = function (params$jscomp$0) {
    params$jscomp$0 = params$jscomp$0 || {};
    this.type = "StateMachine";
    var _super = _jp.Connectors.AbstractBezierConnector.apply(this, arguments);
    var curviness = params$jscomp$0.curviness || 10;
    var margin = params$jscomp$0.margin || 5;
    var proximityLimit = params$jscomp$0.proximityLimit || 80;
    var clockwise =
      params$jscomp$0.orientation &&
      params$jscomp$0.orientation === "clockwise";
    var _controlPoint;
    this._computeBezier = function (paintInfo, params, sp, tp, w, h) {
      var _sx = params.sourcePos[0] < params.targetPos[0] ? 0 : w;
      var _sy = params.sourcePos[1] < params.targetPos[1] ? 0 : h;
      var _tx = params.sourcePos[0] < params.targetPos[0] ? w : 0;
      var _ty = params.sourcePos[1] < params.targetPos[1] ? h : 0;
      if (params.sourcePos[2] === 0) _sx -= margin;
      if (params.sourcePos[2] === 1) _sx += margin;
      if (params.sourcePos[3] === 0) _sy -= margin;
      if (params.sourcePos[3] === 1) _sy += margin;
      if (params.targetPos[2] === 0) _tx -= margin;
      if (params.targetPos[2] === 1) _tx += margin;
      if (params.targetPos[3] === 0) _ty -= margin;
      if (params.targetPos[3] === 1) _ty += margin;
      var _midx = (_sx + _tx) / 2;
      var _midy = (_sy + _ty) / 2;
      var segment = _segment(_sx, _sy, _tx, _ty);
      var distance = Math.sqrt(Math.pow(_tx - _sx, 2) + Math.pow(_ty - _sy, 2));
      _controlPoint = _findControlPoint(
        _midx,
        _midy,
        segment,
        params.sourcePos,
        params.targetPos,
        curviness,
        curviness,
        distance,
        proximityLimit
      );
      var cp1x = _controlPoint[0];
      var cp2x = _controlPoint[0];
      var cp1y = _controlPoint[1];
      var cp2y = _controlPoint[1];
      _super.addSegment(this, "Bezier", {
        x1: _tx,
        y1: _ty,
        x2: _sx,
        y2: _sy,
        cp1x,
        cp1y,
        cp2x,
        cp2y,
      });
    };
  };
  _jp.Connectors.StateMachine = StateMachine;
  _ju.extend(StateMachine, _jp.Connectors.AbstractBezierConnector);
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var STRAIGHT = "Straight";
  var Straight = function (params) {
    this.type = STRAIGHT;
    var _super = _jp.Connectors.AbstractConnector.apply(this, arguments);
    this._compute = function (paintInfo, _) {
      _super.addSegment(this, STRAIGHT, {
        x1: paintInfo.sx,
        y1: paintInfo.sy,
        x2: paintInfo.startStubX,
        y2: paintInfo.startStubY,
      });
      _super.addSegment(this, STRAIGHT, {
        x1: paintInfo.startStubX,
        y1: paintInfo.startStubY,
        x2: paintInfo.endStubX,
        y2: paintInfo.endStubY,
      });
      _super.addSegment(this, STRAIGHT, {
        x1: paintInfo.endStubX,
        y1: paintInfo.endStubY,
        x2: paintInfo.tx,
        y2: paintInfo.ty,
      });
    };
  };
  _jp.Connectors.Straight = Straight;
  _ju.extend(Straight, _jp.Connectors.AbstractConnector);
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var svgAttributeMap = {
    "stroke-linejoin": "stroke-linejoin",
    "stroke-dashoffset": "stroke-dashoffset",
    "stroke-linecap": "stroke-linecap",
  };
  var STROKE_DASHARRAY = "stroke-dasharray";
  var DASHSTYLE = "dashstyle";
  var LINEAR_GRADIENT = "linearGradient";
  var RADIAL_GRADIENT = "radialGradient";
  var DEFS = "defs";
  var FILL = "fill";
  var STOP = "stop";
  var STROKE = "stroke";
  var STROKE_WIDTH = "stroke-width";
  var STYLE = "style";
  var NONE = "none";
  var JSPLUMB_GRADIENT = "jsplumb_gradient_";
  var LINE_WIDTH = "strokeWidth";
  var ns = { svg: "http://www.w3.org/2000/svg" };
  var _attr = function (node, attributes) {
    for (var i in attributes) node.setAttribute(i, "" + attributes[i]);
  };
  var _node = function (name, attributes) {
    attributes = attributes || {};
    attributes.version = "1.1";
    attributes.xmlns = ns.svg;
    return _jp.createElementNS(ns.svg, name, null, null, attributes);
  };
  var _pos = function (d) {
    return "position:absolute;left:" + d[0] + "px;top:" + d[1] + "px";
  };
  var _clearGradient = function (parent) {
    var els = parent.querySelectorAll(" defs,linearGradient,radialGradient");
    for (var i = 0; i < els.length; i++) els[i].parentNode.removeChild(els[i]);
  };
  var _updateGradient = function (
    parent,
    node,
    style,
    dimensions,
    uiComponent
  ) {
    var id = JSPLUMB_GRADIENT + uiComponent._jsPlumb.instance.idstamp();
    _clearGradient(parent);
    if (!style.gradient.offset)
      var g = _node(LINEAR_GRADIENT, { id, gradientUnits: "userSpaceOnUse" });
    else g = _node(RADIAL_GRADIENT, { id });
    var defs = _node(DEFS);
    parent.appendChild(defs);
    defs.appendChild(g);
    for (var i = 0; i < style.gradient.stops.length; i++) {
      var styleToUse =
        uiComponent.segment === 1 || uiComponent.segment === 2
          ? i
          : style.gradient.stops.length - 1 - i;
      var stopColor = style.gradient.stops[styleToUse][1];
      var s = _node(STOP, {
        offset: Math.floor(style.gradient.stops[i][0] * 100) + "%",
        "stop-color": stopColor,
      });
      g.appendChild(s);
    }
    var applyGradientTo = style.stroke ? STROKE : FILL;
    node.setAttribute(applyGradientTo, "url(#" + id + ")");
  };
  var _applyStyles = function (parent, node, style, dimensions, uiComponent) {
    node.setAttribute(FILL, style.fill ? style.fill : NONE);
    node.setAttribute(STROKE, style.stroke ? style.stroke : NONE);
    if (style.gradient)
      _updateGradient(parent, node, style, dimensions, uiComponent);
    else {
      _clearGradient(parent);
      node.setAttribute(STYLE, "");
    }
    if (style.strokeWidth) node.setAttribute(STROKE_WIDTH, style.strokeWidth);
    if (style[DASHSTYLE] && style[LINE_WIDTH] && !style[STROKE_DASHARRAY]) {
      var sep = style[DASHSTYLE].indexOf(",") === -1 ? " " : ",";
      var parts = style[DASHSTYLE].split(sep);
      var styleToUse = "";
      parts.forEach(function (p) {
        styleToUse += Math.floor(p * style.strokeWidth) + sep;
      });
      node.setAttribute(STROKE_DASHARRAY, styleToUse);
    } else if (style[STROKE_DASHARRAY])
      node.setAttribute(STROKE_DASHARRAY, style[STROKE_DASHARRAY]);
    for (var i in svgAttributeMap)
      if (style[i]) node.setAttribute(svgAttributeMap[i], style[i]);
  };
  var _appendAtIndex = function (svg, path, idx) {
    if (svg.childNodes.length > idx)
      svg.insertBefore(path, svg.childNodes[idx]);
    else svg.appendChild(path);
  };
  _ju.svg = { node: _node, attr: _attr, pos: _pos };
  var SvgComponent = function (params) {
    var pointerEventsSpec = params.pointerEventsSpec || "all";
    var renderer = {};
    _jp.jsPlumbUIComponent.apply(this, params.originalArgs);
    this.canvas = null;
    this.path = null;
    this.svg = null;
    this.bgCanvas = null;
    var clazz = params.cssClass + " " + (params.originalArgs[0].cssClass || "");
    var svgParams = {
      style: "",
      width: 0,
      height: 0,
      "pointer-events": pointerEventsSpec,
      position: "absolute",
    };
    this.svg = _node("svg", svgParams);
    if (params.useDivWrapper) {
      this.canvas = _jp.createElement("div", { position: "absolute" });
      _ju.sizeElement(this.canvas, 0, 0, 1, 1);
      this.canvas.className = clazz;
    } else {
      _attr(this.svg, { class: clazz });
      this.canvas = this.svg;
    }
    params._jsPlumb.appendElement(this.canvas, params.originalArgs[0].parent);
    if (params.useDivWrapper) this.canvas.appendChild(this.svg);
    var displayElements = [this.canvas];
    this.getDisplayElements = function () {
      return displayElements;
    };
    this.appendDisplayElement = function (el) {
      displayElements.push(el);
    };
    this.paint = function (style, anchor, extents) {
      if (style != null) {
        var xy = [this.x, this.y];
        var wh = [this.w, this.h];
        if (extents != null) {
          if (extents.xmin < 0) xy[0] += extents.xmin;
          if (extents.ymin < 0) xy[1] += extents.ymin;
          wh[0] = extents.xmax + (extents.xmin < 0 ? -extents.xmin : 0);
          wh[1] = extents.ymax + (extents.ymin < 0 ? -extents.ymin : 0);
        }
        if (params.useDivWrapper) {
          _ju.sizeElement(
            this.canvas,
            xy[0],
            xy[1],
            wh[0] > 0 ? wh[0] : 1,
            wh[1] > 0 ? wh[1] : 1
          );
          xy[0] = 0;
          xy[1] = 0;
          var p = _pos([0, 0]);
        } else p = _pos([xy[0], xy[1]]);
        renderer.paint.apply(this, arguments);
        _attr(this.svg, { style: p, width: wh[0] || 1, height: wh[1] || 1 });
      }
    };
    return { renderer };
  };
  _ju.extend(SvgComponent, _jp.jsPlumbUIComponent, {
    cleanup: function (force) {
      if (force || this.typeId == null) {
        if (this.canvas) this.canvas._jsPlumb = null;
        if (this.svg) this.svg._jsPlumb = null;
        if (this.bgCanvas) this.bgCanvas._jsPlumb = null;
        if (this.canvas && this.canvas.parentNode)
          this.canvas.parentNode.removeChild(this.canvas);
        if (this.bgCanvas && this.bgCanvas.parentNode)
          this.canvas.parentNode.removeChild(this.canvas);
        this.svg = null;
        this.canvas = null;
        this.path = null;
        this.group = null;
      } else {
        if (this.canvas && this.canvas.parentNode)
          this.canvas.parentNode.removeChild(this.canvas);
        if (this.bgCanvas && this.bgCanvas.parentNode)
          this.bgCanvas.parentNode.removeChild(this.bgCanvas);
      }
    },
    reattach: function (instance) {
      var c = instance.getContainer();
      if (this.canvas && this.canvas.parentNode == null)
        c.appendChild(this.canvas);
      if (this.bgCanvas && this.bgCanvas.parentNode == null)
        c.appendChild(this.bgCanvas);
    },
    setVisible: function (v) {
      if (this.canvas) this.canvas.style.display = v ? "block" : "none";
    },
  });
  _jp.ConnectorRenderers.svg = function (params) {
    var self = this;
    var _super = SvgComponent.apply(this, [
      {
        cssClass: params._jsPlumb.connectorClass,
        originalArgs: arguments,
        pointerEventsSpec: "none",
        _jsPlumb: params._jsPlumb,
      },
    ]);
    _super.renderer.paint = function (style, anchor, extents) {
      var segments = self.getSegments();
      var p = "";
      var offset = [0, 0];
      if (extents.xmin < 0) offset[0] = -extents.xmin;
      if (extents.ymin < 0) offset[1] = -extents.ymin;
      if (segments.length > 0) {
        p = self.getPathData();
        var a = {
          d: p,
          transform: "translate(" + offset[0] + "," + offset[1] + ")",
          "pointer-events": params["pointer-events"] || "visibleStroke",
        };
        var outlineStyle = null;
        var d = [self.x, self.y, self.w, self.h];
        if (style.outlineStroke) {
          var outlineWidth = style.outlineWidth || 1;
          var outlineStrokeWidth = style.strokeWidth + 2 * outlineWidth;
          outlineStyle = _jp.extend({}, style);
          delete outlineStyle.gradient;
          outlineStyle.stroke = style.outlineStroke;
          outlineStyle.strokeWidth = outlineStrokeWidth;
          if (self.bgPath == null) {
            self.bgPath = _node("path", a);
            _jp.addClass(self.bgPath, _jp.connectorOutlineClass);
            _appendAtIndex(self.svg, self.bgPath, 0);
          } else _attr(self.bgPath, a);
          _applyStyles(self.svg, self.bgPath, outlineStyle, d, self);
        }
        if (self.path == null) {
          self.path = _node("path", a);
          _appendAtIndex(self.svg, self.path, style.outlineStroke ? 1 : 0);
        } else _attr(self.path, a);
        _applyStyles(self.svg, self.path, style, d, self);
      }
    };
  };
  _ju.extend(_jp.ConnectorRenderers.svg, SvgComponent);
  var SvgEndpoint = (_jp.SvgEndpoint = function (params) {
    var _super = SvgComponent.apply(this, [
      {
        cssClass: params._jsPlumb.endpointClass,
        originalArgs: arguments,
        pointerEventsSpec: "all",
        useDivWrapper: true,
        _jsPlumb: params._jsPlumb,
      },
    ]);
    _super.renderer.paint = function (style) {
      var s = _jp.extend({}, style);
      if (s.outlineStroke) s.stroke = s.outlineStroke;
      if (this.node == null) {
        this.node = this.makeNode(s);
        this.svg.appendChild(this.node);
      } else if (this.updateNode != null) this.updateNode(this.node);
      _applyStyles(
        this.svg,
        this.node,
        s,
        [this.x, this.y, this.w, this.h],
        this
      );
      _pos(this.node, [this.x, this.y]);
    }.bind(this);
  });
  _ju.extend(SvgEndpoint, SvgComponent);
  _jp.Endpoints.svg.Dot = function () {
    _jp.Endpoints.Dot.apply(this, arguments);
    SvgEndpoint.apply(this, arguments);
    this.makeNode = function (style) {
      return _node("circle", {
        cx: this.w / 2,
        cy: this.h / 2,
        r: this.radius,
      });
    };
    this.updateNode = function (node) {
      _attr(node, { cx: this.w / 2, cy: this.h / 2, r: this.radius });
    };
  };
  _ju.extend(_jp.Endpoints.svg.Dot, [_jp.Endpoints.Dot, SvgEndpoint]);
  _jp.Endpoints.svg.Rectangle = function () {
    _jp.Endpoints.Rectangle.apply(this, arguments);
    SvgEndpoint.apply(this, arguments);
    this.makeNode = function (style) {
      return _node("rect", { width: this.w, height: this.h });
    };
    this.updateNode = function (node) {
      _attr(node, { width: this.w, height: this.h });
    };
  };
  _ju.extend(_jp.Endpoints.svg.Rectangle, [
    _jp.Endpoints.Rectangle,
    SvgEndpoint,
  ]);
  _jp.Endpoints.svg.Image = _jp.Endpoints.Image;
  _jp.Endpoints.svg.Blank = _jp.Endpoints.Blank;
  _jp.Overlays.svg.Label = _jp.Overlays.Label;
  _jp.Overlays.svg.Custom = _jp.Overlays.Custom;
  var AbstractSvgArrowOverlay = function (superclass, originalArgs) {
    superclass.apply(this, originalArgs);
    _jp.jsPlumbUIComponent.apply(this, originalArgs);
    this.isAppendedAtTopLevel = false;
    var self = this;
    this.path = null;
    this.paint = function (params, containerExtents) {
      if (params.component.svg && containerExtents) {
        if (this.path == null) {
          this.path = _node("path", { "pointer-events": "all" });
          params.component.svg.appendChild(this.path);
          if (this.elementCreated)
            this.elementCreated(this.path, params.component);
          this.canvas = params.component.svg;
        }
        var clazz =
          originalArgs && originalArgs.length === 1
            ? originalArgs[0].cssClass || ""
            : "";
        var offset = [0, 0];
        if (containerExtents.xmin < 0) offset[0] = -containerExtents.xmin;
        if (containerExtents.ymin < 0) offset[1] = -containerExtents.ymin;
        _attr(this.path, {
          d: makePath(params.d),
          class: clazz,
          stroke: params.stroke ? params.stroke : null,
          fill: params.fill ? params.fill : null,
          transform: "translate(" + offset[0] + "," + offset[1] + ")",
        });
      }
    };
    var makePath = function (d) {
      return isNaN(d.cxy.x) || isNaN(d.cxy.y)
        ? ""
        : "M" +
            d.hxy.x +
            "," +
            d.hxy.y +
            " L" +
            d.tail[0].x +
            "," +
            d.tail[0].y +
            " L" +
            d.cxy.x +
            "," +
            d.cxy.y +
            " L" +
            d.tail[1].x +
            "," +
            d.tail[1].y +
            " L" +
            d.hxy.x +
            "," +
            d.hxy.y;
    };
    this.transfer = function (target) {
      if (target.canvas && this.path && this.path.parentNode) {
        this.path.parentNode.removeChild(this.path);
        target.canvas.appendChild(this.path);
      }
    };
  };
  var svgProtoFunctions = {
    cleanup: function (force) {
      if (this.path != null)
        if (force) this._jsPlumb.instance.removeElement(this.path);
        else if (this.path.parentNode)
          this.path.parentNode.removeChild(this.path);
    },
    reattach: function (instance, component) {
      if (this.path && component.canvas)
        component.canvas.appendChild(this.path);
    },
    setVisible: function (v) {
      if (this.path != null) this.path.style.display = v ? "block" : "none";
    },
  };
  _ju.extend(AbstractSvgArrowOverlay, [
    _jp.jsPlumbUIComponent,
    _jp.Overlays.AbstractOverlay,
  ]);
  _jp.Overlays.svg.Arrow = function () {
    AbstractSvgArrowOverlay.apply(this, [_jp.Overlays.Arrow, arguments]);
  };
  _ju.extend(
    _jp.Overlays.svg.Arrow,
    [_jp.Overlays.Arrow, AbstractSvgArrowOverlay],
    svgProtoFunctions
  );
  _jp.Overlays.svg.PlainArrow = function () {
    AbstractSvgArrowOverlay.apply(this, [_jp.Overlays.PlainArrow, arguments]);
  };
  _ju.extend(
    _jp.Overlays.svg.PlainArrow,
    [_jp.Overlays.PlainArrow, AbstractSvgArrowOverlay],
    svgProtoFunctions
  );
  _jp.Overlays.svg.Diamond = function () {
    AbstractSvgArrowOverlay.apply(this, [_jp.Overlays.Diamond, arguments]);
  };
  _ju.extend(
    _jp.Overlays.svg.Diamond,
    [_jp.Overlays.Diamond, AbstractSvgArrowOverlay],
    svgProtoFunctions
  );
  _jp.Overlays.svg.GuideLines = function () {
    var path = null;
    var self = this;
    var p1_1;
    var p1_2;
    _jp.Overlays.GuideLines.apply(this, arguments);
    this.paint = function (params, containerExtents) {
      if (path == null) {
        path = _node("path");
        params.connector.svg.appendChild(path);
        self.attachListeners(path, params.connector);
        self.attachListeners(path, self);
        p1_1 = _node("path");
        params.connector.svg.appendChild(p1_1);
        self.attachListeners(p1_1, params.connector);
        self.attachListeners(p1_1, self);
        p1_2 = _node("path");
        params.connector.svg.appendChild(p1_2);
        self.attachListeners(p1_2, params.connector);
        self.attachListeners(p1_2, self);
      }
      var offset = [0, 0];
      if (containerExtents.xmin < 0) offset[0] = -containerExtents.xmin;
      if (containerExtents.ymin < 0) offset[1] = -containerExtents.ymin;
      _attr(path, {
        d: makePath(params.head, params.tail),
        stroke: "red",
        fill: null,
        transform: "translate(" + offset[0] + "," + offset[1] + ")",
      });
      _attr(p1_1, {
        d: makePath(params.tailLine[0], params.tailLine[1]),
        stroke: "blue",
        fill: null,
        transform: "translate(" + offset[0] + "," + offset[1] + ")",
      });
      _attr(p1_2, {
        d: makePath(params.headLine[0], params.headLine[1]),
        stroke: "green",
        fill: null,
        transform: "translate(" + offset[0] + "," + offset[1] + ")",
      });
    };
    var makePath = function (d1, d2) {
      return "M " + d1.x + "," + d1.y + " L" + d2.x + "," + d2.y;
    };
  };
  _ju.extend(_jp.Overlays.svg.GuideLines, _jp.Overlays.GuideLines);
}).call(typeof window !== "undefined" ? window : this);
(function () {
  var root = this;
  var _jp = root.jsPlumb;
  var _ju = root.jsPlumbUtil;
  var _jk = root.Katavorio;
  var _jg = root.Biltong;
  var _getEventManager = function (instance) {
    var e = instance._mottle;
    if (!e) e = instance._mottle = new root.Mottle();
    return e;
  };
  var _getDragManager = function (instance, category) {
    category = category || "main";
    var key = "_katavorio_" + category;
    var k = instance[key];
    var e = instance.getEventManager();
    if (!k) {
      k = new _jk({
        bind: e.on,
        unbind: e.off,
        getSize: _jp.getSize,
        getConstrainingRectangle: function (el) {
          return [el.parentNode.scrollWidth, el.parentNode.scrollHeight];
        },
        getPosition: function (el, relativeToRoot) {
          var o = instance.getOffset(
            el,
            relativeToRoot,
            el._katavorioDrag ? el.offsetParent : null
          );
          return [o.left, o.top];
        },
        setPosition: function (el, xy) {
          el.style.left = xy[0] + "px";
          el.style.top = xy[1] + "px";
        },
        addClass: _jp.addClass,
        removeClass: _jp.removeClass,
        intersects: _jg.intersects,
        indexOf: function (l, i) {
          return l.indexOf(i);
        },
        scope: instance.getDefaultScope(),
        css: {
          noSelect: instance.dragSelectClass,
          droppable: "jtk-droppable",
          draggable: "jtk-draggable",
          drag: "jtk-drag",
          selected: "jtk-drag-selected",
          active: "jtk-drag-active",
          hover: "jtk-drag-hover",
          ghostProxy: "jtk-ghost-proxy",
        },
      });
      k.setZoom(instance.getZoom());
      instance[key] = k;
      instance.bind("zoom", k.setZoom);
    }
    return k;
  };
  var _dragStart = function (params) {
    var options = params.el._jsPlumbDragOptions;
    var cont = true;
    if (options.canDrag) cont = options.canDrag();
    if (cont) {
      this.setHoverSuspended(true);
      this.select({ source: params.el }).addClass(
        this.elementDraggingClass + " " + this.sourceElementDraggingClass,
        true
      );
      this.select({ target: params.el }).addClass(
        this.elementDraggingClass + " " + this.targetElementDraggingClass,
        true
      );
      this.setConnectionBeingDragged(true);
    }
    return cont;
  };
  var _dragMove = function (params) {
    var ui = this.getUIPosition(arguments, this.getZoom());
    if (ui != null) {
      var o = params.el._jsPlumbDragOptions;
      this.draw(params.el, ui, null, true);
      if (o._dragging) this.addClass(params.el, "jtk-dragged");
      o._dragging = true;
    }
  };
  var _dragStop = function (params) {
    var elements = params.selection;
    var uip;
    var _one = function (_e) {
      if (_e[1] != null) {
        uip = this.getUIPosition([
          { el: _e[2].el, pos: [_e[1].left, _e[1].top] },
        ]);
        var drawResult = this.draw(_e[2].el, uip);
      }
      if (_e[0]._jsPlumbDragOptions != null)
        delete _e[0]._jsPlumbDragOptions._dragging;
      this.removeClass(_e[0], "jtk-dragged");
      this.select({ source: _e[2].el }).removeClass(
        this.elementDraggingClass + " " + this.sourceElementDraggingClass,
        true
      );
      this.select({ target: _e[2].el }).removeClass(
        this.elementDraggingClass + " " + this.targetElementDraggingClass,
        true
      );
      params.e._drawResult = params.e._drawResult || { c: [], e: [], a: [] };
      Array.prototype.push.apply(params.e._drawResult.c, drawResult.c);
      Array.prototype.push.apply(params.e._drawResult.e, drawResult.e);
      Array.prototype.push.apply(params.e._drawResult.a, drawResult.a);
      this.getDragManager().dragEnded(_e[2].el);
    }.bind(this);
    for (var i = 0; i < elements.length; i++) _one(elements[i]);
    this.setHoverSuspended(false);
    this.setConnectionBeingDragged(false);
  };
  var _animProps = function (o, p) {
    var _one = function (pName) {
      if (p[pName] != null)
        if (_ju.isString(p[pName])) {
          var m = p[pName].match(/-=/) ? -1 : 1;
          var v = p[pName].substring(2);
          return o[pName] + m * v;
        } else return p[pName];
      else return o[pName];
    };
    return [_one("left"), _one("top")];
  };
  var _genLoc = function (prefix, e) {
    if (e == null) return [0, 0];
    var ts = _touches(e);
    var t = _getTouch(ts, 0);
    return [t[prefix + "X"], t[prefix + "Y"]];
  };
  var _pageLocation = _genLoc.bind(this, "page");
  var _screenLocation = _genLoc.bind(this, "screen");
  var _clientLocation = _genLoc.bind(this, "client");
  var _getTouch = function (touches, idx) {
    return touches.item ? touches.item(idx) : touches[idx];
  };
  var _touches = function (e) {
    return e.touches && e.touches.length > 0
      ? e.touches
      : e.changedTouches && e.changedTouches.length > 0
      ? e.changedTouches
      : e.targetTouches && e.targetTouches.length > 0
      ? e.targetTouches
      : [e];
  };
  var DragManager = function (_currentInstance) {
    var _draggables = {};
    var _dlist = [];
    var _delements = {};
    var _elementsWithEndpoints = {};
    var _draggablesForElements = {};
    this.register = function (el) {
      var id = _currentInstance.getId(el);
      var parentOffset;
      if (!_draggables[id]) {
        _draggables[id] = el;
        _dlist.push(el);
        _delements[id] = {};
      }
      var _oneLevel = function (p) {
        if (p)
          for (var i = 0; i < p.childNodes.length; i++)
            if (
              p.childNodes[i].nodeType !== 3 &&
              p.childNodes[i].nodeType !== 8
            ) {
              var cEl = jsPlumb.getElement(p.childNodes[i]);
              var cid = _currentInstance.getId(p.childNodes[i], null, true);
              if (
                cid &&
                _elementsWithEndpoints[cid] &&
                _elementsWithEndpoints[cid] > 0
              ) {
                if (!parentOffset)
                  parentOffset = _currentInstance.getOffset(el);
                var cOff = _currentInstance.getOffset(cEl);
                _delements[id][cid] = {
                  id: cid,
                  offset: {
                    left: cOff.left - parentOffset.left,
                    top: cOff.top - parentOffset.top,
                  },
                };
                _draggablesForElements[cid] = id;
              }
              _oneLevel(p.childNodes[i]);
            }
      };
      _oneLevel(el);
    };
    this.updateOffsets = function (elId, childOffsetOverrides) {
      if (elId != null) {
        childOffsetOverrides = childOffsetOverrides || {};
        var domEl = jsPlumb.getElement(elId);
        var id = _currentInstance.getId(domEl);
        var children = _delements[id];
        if (children)
          for (var i in children)
            if (children.hasOwnProperty(i)) {
              var cel = jsPlumb.getElement(i);
              var cOff =
                childOffsetOverrides[i] || _currentInstance.getOffset(cel);
              if (cel.offsetParent == null && _delements[id][i] != null)
                continue;
              if (!parentOffset)
                var parentOffset = _currentInstance.getOffset(domEl);
              _delements[id][i] = {
                id: i,
                offset: {
                  left: cOff.left - parentOffset.left,
                  top: cOff.top - parentOffset.top,
                },
              };
              _draggablesForElements[i] = id;
            }
      }
    };
    this.endpointAdded = function (el, id) {
      id = id || _currentInstance.getId(el);
      var b = document.body;
      var p = el.parentNode;
      for (
        _elementsWithEndpoints[id] = _elementsWithEndpoints[id]
          ? _elementsWithEndpoints[id] + 1
          : 1;
        p != null && p !== b;

      ) {
        var pid = _currentInstance.getId(p, null, true);
        if (pid && _draggables[pid]) {
          var pLoc = _currentInstance.getOffset(p);
          if (_delements[pid][id] == null) {
            var cLoc = _currentInstance.getOffset(el);
            _delements[pid][id] = {
              id,
              offset: { left: cLoc.left - pLoc.left, top: cLoc.top - pLoc.top },
            };
            _draggablesForElements[id] = pid;
          }
          break;
        }
        p = p.parentNode;
      }
    };
    this.endpointDeleted = function (endpoint) {
      if (_elementsWithEndpoints[endpoint.elementId]) {
        _elementsWithEndpoints[endpoint.elementId]--;
        if (_elementsWithEndpoints[endpoint.elementId] <= 0)
          for (var i in _delements)
            if (_delements.hasOwnProperty(i) && _delements[i]) {
              delete _delements[i][endpoint.elementId];
              delete _draggablesForElements[endpoint.elementId];
            }
      }
    };
    this.changeId = function (oldId, newId) {
      _delements[newId] = _delements[oldId];
      _delements[oldId] = {};
      _draggablesForElements[newId] = _draggablesForElements[oldId];
      _draggablesForElements[oldId] = null;
    };
    this.getElementsForDraggable = function (id) {
      return _delements[id];
    };
    this.elementRemoved = function (elementId) {
      var elId = _draggablesForElements[elementId];
      if (elId) {
        _delements[elId] && delete _delements[elId][elementId];
        delete _draggablesForElements[elementId];
      }
    };
    this.reset = function () {
      _draggables = {};
      _dlist = [];
      _delements = {};
      _elementsWithEndpoints = {};
    };
    this.dragEnded = function (el) {
      if (el.offsetParent != null) {
        var id = _currentInstance.getId(el);
        var ancestor = _draggablesForElements[id];
        if (ancestor) this.updateOffsets(ancestor);
      }
    };
    this.setParent = function (el, elId, p, pId, currentChildLocation) {
      var current = _draggablesForElements[elId];
      if (!_delements[pId]) _delements[pId] = {};
      var pLoc = _currentInstance.getOffset(p);
      var cLoc = currentChildLocation || _currentInstance.getOffset(el);
      if (current && _delements[current]) delete _delements[current][elId];
      _delements[pId][elId] = {
        id: elId,
        offset: { left: cLoc.left - pLoc.left, top: cLoc.top - pLoc.top },
      };
      _draggablesForElements[elId] = pId;
    };
    this.clearParent = function (el, elId) {
      var current = _draggablesForElements[elId];
      if (current) {
        delete _delements[current][elId];
        delete _draggablesForElements[elId];
      }
    };
    this.revalidateParent = function (el, elId, childOffset) {
      var current = _draggablesForElements[elId];
      if (current) {
        var co = {};
        co[elId] = childOffset;
        this.updateOffsets(current, co);
        _currentInstance.revalidate(current);
      }
    };
    this.getDragAncestor = function (el) {
      var de = jsPlumb.getElement(el);
      var id = _currentInstance.getId(de);
      var aid = _draggablesForElements[id];
      if (aid) return jsPlumb.getElement(aid);
      else return null;
    };
  };
  var _setClassName = function (el, cn, classList) {
    cn = _ju.fastTrim(cn);
    if (typeof el.className.baseVal !== "undefined") el.className.baseVal = cn;
    else el.className = cn;
    try {
      var cl = el.classList;
      if (cl != null) {
        for (; cl.length > 0; ) cl.remove(cl.item(0));
        for (var i = 0; i < classList.length; i++)
          if (classList[i]) cl.add(classList[i]);
      }
    } catch (e) {
      _ju.log("JSPLUMB: cannot set class list", e);
    }
  };
  var _getClassName = function (el) {
    return typeof el.className.baseVal === "undefined"
      ? el.className
      : el.className.baseVal;
  };
  var _classManip = function (el, classesToAdd, classesToRemove) {
    classesToAdd =
      classesToAdd == null
        ? []
        : _ju.isArray(classesToAdd)
        ? classesToAdd
        : classesToAdd.split(/\s+/);
    classesToRemove =
      classesToRemove == null
        ? []
        : _ju.isArray(classesToRemove)
        ? classesToRemove
        : classesToRemove.split(/\s+/);
    var className = _getClassName(el);
    var curClasses = className.split(/\s+/);
    var _oneSet = function (add, classes) {
      for (var i = 0; i < classes.length; i++)
        if (add) {
          if (curClasses.indexOf(classes[i]) === -1)
            curClasses.push(classes[i]);
        } else {
          var idx = curClasses.indexOf(classes[i]);
          if (idx !== -1) curClasses.splice(idx, 1);
        }
    };
    _oneSet(true, classesToAdd);
    _oneSet(false, classesToRemove);
    _setClassName(el, curClasses.join(" "), curClasses);
  };
  root.jsPlumb.extend(root.jsPlumbInstance.prototype, {
    headless: false,
    pageLocation: _pageLocation,
    screenLocation: _screenLocation,
    clientLocation: _clientLocation,
    getDragManager: function () {
      if (this.dragManager == null) this.dragManager = new DragManager(this);
      return this.dragManager;
    },
    recalculateOffsets: function (elId) {
      this.getDragManager().updateOffsets(elId);
    },
    createElement: function (tag, style, clazz, atts) {
      return this.createElementNS(null, tag, style, clazz, atts);
    },
    createElementNS: function (ns, tag, style, clazz, atts) {
      var e =
        ns == null
          ? document.createElement(tag)
          : document.createElementNS(ns, tag);
      var i;
      style = style || {};
      for (i in style) e.style[i] = style[i];
      if (clazz) e.className = clazz;
      atts = atts || {};
      for (i in atts) e.setAttribute(i, "" + atts[i]);
      return e;
    },
    getAttribute: function (el, attName) {
      return el.getAttribute != null ? el.getAttribute(attName) : null;
    },
    setAttribute: function (el, a, v) {
      if (el.setAttribute != null) el.setAttribute(a, v);
    },
    setAttributes: function (el, atts) {
      for (var i in atts)
        if (atts.hasOwnProperty(i)) el.setAttribute(i, atts[i]);
    },
    appendToRoot: function (node) {
      document.body.appendChild(node);
    },
    getRenderModes: function () {
      return ["svg"];
    },
    getClass: _getClassName,
    addClass: function (el, clazz) {
      jsPlumb.each(el, function (e) {
        _classManip(e, clazz);
      });
    },
    hasClass: function (el, clazz) {
      el = jsPlumb.getElement(el);
      if (el.classList) return el.classList.contains(clazz);
      else return _getClassName(el).indexOf(clazz) !== -1;
    },
    removeClass: function (el, clazz) {
      jsPlumb.each(el, function (e) {
        _classManip(e, null, clazz);
      });
    },
    toggleClass: function (el, clazz) {
      if (jsPlumb.hasClass(el, clazz)) jsPlumb.removeClass(el, clazz);
      else jsPlumb.addClass(el, clazz);
    },
    updateClasses: function (el, toAdd, toRemove) {
      jsPlumb.each(el, function (e) {
        _classManip(e, toAdd, toRemove);
      });
    },
    setClass: function (el, clazz) {
      if (clazz != null)
        jsPlumb.each(el, function (e) {
          _setClassName(e, clazz, clazz.split(/\s+/));
        });
    },
    setPosition: function (el, p) {
      el.style.left = p.left + "px";
      el.style.top = p.top + "px";
    },
    getPosition: function (el) {
      var _one = function (prop) {
        var v = el.style[prop];
        return v ? v.substring(0, v.length - 2) : 0;
      };
      return { left: _one("left"), top: _one("top") };
    },
    getStyle: function (el, prop) {
      if (typeof window.getComputedStyle !== "undefined")
        return getComputedStyle(el, null).getPropertyValue(prop);
      else return el.currentStyle[prop];
    },
    getSelector: function (ctx, spec) {
      var sel = null;
      if (arguments.length === 1)
        sel = ctx.nodeType != null ? ctx : document.querySelectorAll(ctx);
      else sel = ctx.querySelectorAll(spec);
      return sel;
    },
    getOffset: function (el, relativeToRoot, container) {
      el = jsPlumb.getElement(el);
      container = container || this.getContainer();
      var out = { left: el.offsetLeft, top: el.offsetTop };
      var op =
        relativeToRoot ||
        (container != null && el !== container && el.offsetParent !== container)
          ? el.offsetParent
          : null;
      for (
        var _maybeAdjustScroll = function (offsetParent) {
          if (
            offsetParent != null &&
            offsetParent !== document.body &&
            (offsetParent.scrollTop > 0 || offsetParent.scrollLeft > 0)
          ) {
            out.left -= offsetParent.scrollLeft;
            out.top -= offsetParent.scrollTop;
          }
        }.bind(this);
        op != null;

      ) {
        out.left += op.offsetLeft;
        out.top += op.offsetTop;
        _maybeAdjustScroll(op);
        op = relativeToRoot
          ? op.offsetParent
          : op.offsetParent === container
          ? null
          : op.offsetParent;
      }
      if (
        container != null &&
        !relativeToRoot &&
        (container.scrollTop > 0 || container.scrollLeft > 0)
      ) {
        var pp =
          el.offsetParent != null
            ? this.getStyle(el.offsetParent, "position")
            : "static";
        var p = this.getStyle(el, "position");
        if (
          p !== "absolute" &&
          p !== "fixed" &&
          pp !== "absolute" &&
          pp !== "fixed"
        ) {
          out.left -= container.scrollLeft;
          out.top -= container.scrollTop;
        }
      }
      return out;
    },
    getPositionOnElement: function (evt, el, zoom) {
      var box =
        typeof el.getBoundingClientRect !== "undefined"
          ? el.getBoundingClientRect()
          : { left: 0, top: 0, width: 0, height: 0 };
      var body = document.body;
      var docElem = document.documentElement;
      var scrollTop = window.pageYOffset || docElem.scrollTop || body.scrollTop;
      var scrollLeft =
        window.pageXOffset || docElem.scrollLeft || body.scrollLeft;
      var clientTop = docElem.clientTop || body.clientTop || 0;
      var clientLeft = docElem.clientLeft || body.clientLeft || 0;
      var pst = 0;
      var psl = 0;
      var top = box.top + scrollTop - clientTop + pst * zoom;
      var left = box.left + scrollLeft - clientLeft + psl * zoom;
      var cl = jsPlumb.pageLocation(evt);
      var w = box.width || el.offsetWidth * zoom;
      var h = box.height || el.offsetHeight * zoom;
      var x = (cl[0] - left) / w;
      var y = (cl[1] - top) / h;
      return [x, y];
    },
    getAbsolutePosition: function (el) {
      var _one = function (s) {
        var ss = el.style[s];
        if (ss) return parseFloat(ss.substring(0, ss.length - 2));
      };
      return [_one("left"), _one("top")];
    },
    setAbsolutePosition: function (el, xy, animateFrom, animateOptions) {
      if (animateFrom)
        this.animate(
          el,
          {
            left: "+\x3d" + (xy[0] - animateFrom[0]),
            top: "+\x3d" + (xy[1] - animateFrom[1]),
          },
          animateOptions
        );
      else {
        el.style.left = xy[0] + "px";
        el.style.top = xy[1] + "px";
      }
    },
    getSize: function (el) {
      return [el.offsetWidth, el.offsetHeight];
    },
    getWidth: function (el) {
      return el.offsetWidth;
    },
    getHeight: function (el) {
      return el.offsetHeight;
    },
    getRenderMode: function () {
      return "svg";
    },
    draggable: function (el, options) {
      var info;
      el =
        _ju.isArray(el) || (el.length != null && !_ju.isString(el)) ? el : [el];
      Array.prototype.slice.call(el).forEach(
        function (_el) {
          info = this.info(_el);
          if (info.el)
            this._initDraggableIfNecessary(
              info.el,
              true,
              options,
              info.id,
              true
            );
        }.bind(this)
      );
      return this;
    },
    snapToGrid: function (el, x, y) {
      var out = [];
      var _oneEl = function (_el) {
        var info = this.info(_el);
        if (info.el != null && info.el._katavorioDrag) {
          var snapped = info.el._katavorioDrag.snap(x, y);
          this.revalidate(info.el);
          out.push([info.el, snapped]);
        }
      }.bind(this);
      if (arguments.length === 1 || arguments.length === 3) _oneEl(el, x, y);
      else {
        var _me = this.getManagedElements();
        for (var mel in _me) _oneEl(mel, arguments[0], arguments[1]);
      }
      return out;
    },
    initDraggable: function (el, options, category) {
      _getDragManager(this, category).draggable(el, options);
      el._jsPlumbDragOptions = options;
    },
    destroyDraggable: function (el, category) {
      _getDragManager(this, category).destroyDraggable(el);
      delete el._jsPlumbDragOptions;
    },
    unbindDraggable: function (el, evt, fn, category) {
      _getDragManager(this, category).destroyDraggable(el, evt, fn);
    },
    setDraggable: function (element, draggable) {
      return jsPlumb.each(
        element,
        function (el) {
          if (this.isDragSupported(el)) {
            this._draggableStates[this.getAttribute(el, "id")] = draggable;
            this.setElementDraggable(el, draggable);
          }
        }.bind(this)
      );
    },
    _draggableStates: {},
    toggleDraggable: function (el$jscomp$0) {
      var state;
      jsPlumb.each(
        el$jscomp$0,
        function (el) {
          var elId = this.getAttribute(el, "id");
          state =
            this._draggableStates[elId] == null
              ? false
              : this._draggableStates[elId];
          state = !state;
          this._draggableStates[elId] = state;
          this.setDraggable(el, state);
          return state;
        }.bind(this)
      );
      return state;
    },
    _initDraggableIfNecessary: function (
      element,
      isDraggable,
      dragOptions,
      id,
      fireEvent
    ) {
      if (!jsPlumb.headless) {
        var _draggable = isDraggable == null ? false : isDraggable;
        if (_draggable)
          if (jsPlumb.isDragSupported(element, this)) {
            var options = dragOptions || this.Defaults.DragOptions;
            options = jsPlumb.extend({}, options);
            if (!jsPlumb.isAlreadyDraggable(element, this)) {
              var dragEvent = jsPlumb.dragEvents.drag;
              var stopEvent = jsPlumb.dragEvents.stop;
              var startEvent = jsPlumb.dragEvents.start;
              this.manage(id, element);
              options[startEvent] = _ju.wrap(
                options[startEvent],
                _dragStart.bind(this)
              );
              options[dragEvent] = _ju.wrap(
                options[dragEvent],
                _dragMove.bind(this)
              );
              options[stopEvent] = _ju.wrap(
                options[stopEvent],
                _dragStop.bind(this)
              );
              var elId = this.getId(element);
              this._draggableStates[elId] = true;
              var draggable = this._draggableStates[elId];
              options.disabled = draggable == null ? false : !draggable;
              this.initDraggable(element, options);
              this.getDragManager().register(element);
              if (fireEvent)
                this.fire("elementDraggable", { el: element, options });
            } else if (dragOptions.force) this.initDraggable(element, options);
          }
      }
    },
    animationSupported: true,
    getElement: function (el) {
      if (el == null) return null;
      el =
        typeof el === "string"
          ? el
          : el.tagName == null && el.length != null && el.enctype == null
          ? el[0]
          : el;
      return typeof el === "string" ? document.getElementById(el) : el;
    },
    removeElement: function (element) {
      _getDragManager(this).elementRemoved(element);
      this.getEventManager().remove(element);
    },
    doAnimate: function (el, properties, options) {
      options = options || {};
      var o = this.getOffset(el);
      var ap = _animProps(o, properties);
      var ldist = ap[0] - o.left;
      var tdist = ap[1] - o.top;
      var d = options.duration || 250;
      var step = 15;
      var steps = d / step;
      var linc = (step / d) * ldist;
      var tinc = (step / d) * tdist;
      var idx = 0;
      var _int = setInterval(function () {
        _jp.setPosition(el, {
          left: o.left + linc * (idx + 1),
          top: o.top + tinc * (idx + 1),
        });
        if (options.step != null) options.step(idx, Math.ceil(steps));
        idx++;
        if (idx >= steps) {
          window.clearInterval(_int);
          if (options.complete != null) options.complete();
        }
      }, step);
    },
    destroyDroppable: function (el, category) {
      _getDragManager(this, category).destroyDroppable(el);
    },
    unbindDroppable: function (el, evt, fn, category) {
      _getDragManager(this, category).destroyDroppable(el, evt, fn);
    },
    droppable: function (el, options) {
      el =
        _ju.isArray(el) || (el.length != null && !_ju.isString(el)) ? el : [el];
      var info;
      options = options || {};
      options.allowLoopback = false;
      Array.prototype.slice.call(el).forEach(
        function (_el) {
          info = this.info(_el);
          if (info.el) this.initDroppable(info.el, options);
        }.bind(this)
      );
      return this;
    },
    initDroppable: function (el, options, category) {
      _getDragManager(this, category).droppable(el, options);
    },
    isAlreadyDraggable: function (el) {
      return el._katavorioDrag != null;
    },
    isDragSupported: function (el, options) {
      return true;
    },
    isDropSupported: function (el, options) {
      return true;
    },
    isElementDraggable: function (el) {
      el = _jp.getElement(el);
      return el._katavorioDrag && el._katavorioDrag.isEnabled();
    },
    getDragObject: function (eventArgs) {
      return eventArgs[0].drag.getDragElement();
    },
    getDragScope: function (el) {
      return (el._katavorioDrag && el._katavorioDrag.scopes.join(" ")) || "";
    },
    getDropEvent: function (args) {
      return args[0].e;
    },
    getUIPosition: function (eventArgs, zoom) {
      var el = eventArgs[0].el;
      if (el.offsetParent == null) return null;
      var finalPos = eventArgs[0].finalPos || eventArgs[0].pos;
      var p = { left: finalPos[0], top: finalPos[1] };
      if (el._katavorioDrag && el.offsetParent !== this.getContainer()) {
        var oc = this.getOffset(el.offsetParent);
        p.left += oc.left;
        p.top += oc.top;
      }
      return p;
    },
    setDragFilter: function (el, filter, _exclude) {
      if (el._katavorioDrag) el._katavorioDrag.setFilter(filter, _exclude);
    },
    setElementDraggable: function (el, draggable) {
      el = _jp.getElement(el);
      if (el._katavorioDrag) el._katavorioDrag.setEnabled(draggable);
    },
    setDragScope: function (el, scope) {
      if (el._katavorioDrag) el._katavorioDrag.k.setDragScope(el, scope);
    },
    setDropScope: function (el, scope) {
      if (el._katavorioDrop && el._katavorioDrop.length > 0)
        el._katavorioDrop[0].k.setDropScope(el, scope);
    },
    addToPosse: function (el, spec) {
      var specs = Array.prototype.slice.call(arguments, 1);
      var dm = _getDragManager(this);
      _jp.each(el, function (_el) {
        _el = [_jp.getElement(_el)];
        _el.push.apply(_el, specs);
        dm.addToPosse.apply(dm, _el);
      });
    },
    setPosse: function (el, spec) {
      var specs = Array.prototype.slice.call(arguments, 1);
      var dm = _getDragManager(this);
      _jp.each(el, function (_el) {
        _el = [_jp.getElement(_el)];
        _el.push.apply(_el, specs);
        dm.setPosse.apply(dm, _el);
      });
    },
    removeFromPosse: function (el, posseId) {
      var specs = Array.prototype.slice.call(arguments, 1);
      var dm = _getDragManager(this);
      _jp.each(el, function (_el) {
        _el = [_jp.getElement(_el)];
        _el.push.apply(_el, specs);
        dm.removeFromPosse.apply(dm, _el);
      });
    },
    removeFromAllPosses: function (el) {
      var dm = _getDragManager(this);
      _jp.each(el, function (_el) {
        dm.removeFromAllPosses(_jp.getElement(_el));
      });
    },
    setPosseState: function (el, posseId, state) {
      var dm = _getDragManager(this);
      _jp.each(el, function (_el) {
        dm.setPosseState(_jp.getElement(_el), posseId, state);
      });
    },
    dragEvents: {
      start: "start",
      stop: "stop",
      drag: "drag",
      step: "step",
      over: "over",
      out: "out",
      drop: "drop",
      complete: "complete",
      beforeStart: "beforeStart",
    },
    animEvents: { step: "step", complete: "complete" },
    stopDrag: function (el) {
      if (el._katavorioDrag) el._katavorioDrag.abort();
    },
    addToDragSelection: function (spec) {
      var el = this.getElement(spec);
      if (el != null && (el._isJsPlumbGroup || el._jsPlumbGroup == null))
        _getDragManager(this).select(spec);
    },
    removeFromDragSelection: function (spec) {
      _getDragManager(this).deselect(spec);
    },
    getDragSelection: function () {
      return _getDragManager(this).getSelection();
    },
    clearDragSelection: function () {
      _getDragManager(this).deselectAll();
    },
    trigger: function (el, event, originalEvent, payload) {
      this.getEventManager().trigger(el, event, originalEvent, payload);
    },
    doReset: function () {
      for (var key in this)
        if (key.indexOf("_katavorio_") === 0) this[key].reset();
    },
    getEventManager: function () {
      return _getEventManager(this);
    },
    on: function (el, event, callback) {
      this.getEventManager().on.apply(this, arguments);
      return this;
    },
    off: function (el, event, callback) {
      this.getEventManager().off.apply(this, arguments);
      return this;
    },
  });
  var ready = function (f) {
    var _do = function () {
      if (
        /complete|loaded|interactive/.test(document.readyState) &&
        typeof document.body !== "undefined" &&
        document.body != null
      )
        f();
      else setTimeout(_do, 9);
    };
    _do();
  };
  ready(_jp.init);
}).call(typeof window !== "undefined" ? window : this);

// 2.15.6
