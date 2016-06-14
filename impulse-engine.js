class Vec2 {
  constructor(x, y) {
    this.set(x, y);
  }
  set(x = 0, y = 0) {
    if (typeof x.x === 'number') {
      this.x = x.x;
      this.y = x.y;
    } else {
      this.x = x;
      this.y = y;
    }
    return this;
  }
  clone() {
    return new Vec2(this);
  }
  negi() {
    return this.neg(this);
  }
  neg(out = new Vec2()) {
    out.x = -this.x;
    out.y = -this.y;
    return out;
  }
  muli(s) {
    return this.mul(s, this);
  }
  mul(s, out = new Vec2()) {
    if (typeof s.x !== 'number') {
      out.x = s * this.x;
      out.y = s * this.y;
    } else {
      out.x = s.x * this.x;
      out.y = s.y * this.y;
    }
    return out;
  }
  divi(s) {
    return this.div(s, this);
  }
  div(s, out = new Vec2()) {
    if (typeof s.x !== 'number') {
      out.x = this.x / s;
      out.y = this.y / s;
    } else {
      out.x = this.x / s.x;
      out.y = this.y / s.y;
    }
    return out;
  }
  addi(s) {
    return this.add(s, this);
  }
  add(s, out = new Vec2()) {
    if (typeof s.x !== 'number') {
      out.x = this.x + s;
      out.y = this.y + s;
    } else {
      out.x = this.x + s.x;
      out.y = this.y + s.y;
    }
    return this;
  }
  addsi(v, s) {
    return this.adds(v, s, this);
  }
  adds(v, s, out = new Vec2()) {
    out.x = this.x + v.x * s;
    out.y = this.y + v.y * s;
    return out;
  }
  subi(v) {
    return this.sub(v, this);
  }
  sub(v, out = new Vec2()) {
    if (typeof v.x !== 'number') {
      out.x = this.x - v;
      out.y = this.y - v;
    } else {
      out.x = this.x - v.x;
      out.y = this.y - v.y;
    }
    return out;
  }
  lengthSq() {
    return this.x * this.x + this.y * this.y;
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y);
  }
  rotate(radians) {
    const c = Math.cos(radians);
    const s = Math.sin(radians);
    this.x = this.x * c - this.y * s;
    this.y = this.x * s + this.y * c;
    return this;
  }
  normalize() {
    const lenSq = this.lengthSq();
    if (lenSq > ImpulseMath.EPSILON_SQ) {
      const invLen = 1 / Math.sqrt(lenSq);
      this.x *= invLen;
      this.y *= invLen;
    }
    return this;
  }
  mini(a, b) {
    return this.min(a, b, this);
  }
  maxi(a, b) {
    return this.max(a, b, this);
  }
  min(a, b, out = new Vec2()) {
    out.x = Math.min(a.x, b.x);
    out.y = Math.min(a.y, b.y);
    return out;
  }
  max(a, b, out = new Vec2()) {
    out.x = Math.max(a.x, b.x);
    out.y = Math.max(a.y, b.y);
    return out;
  }
  dot(v) {
    return this.x * v.x + this.y * v.y;
  }
  distanceSq(v) {
    const dx = this.x - v.x;
    const dy = this.y - v.y;
    return dx * dx + dy * dy;
  }
  distance(v) {
    const dx = this.x - v.x;
    const dy = this.y - v.y;
    return Math.sqrt(dx * dx + dy * dy);
  }
  crossVS(v, a) {
    this.x = v.y * a;
    this.y = v.x * -a;
    return this;
  }
  crossSV(a, v) {
    this.x = v.y * -a;
    this.y = v.x * a;
    return this;
  }
  crossV(v) {
    return this.x * v.y - this.y * v.x;
  }
  static arrayOf(length) {
    const array = new Array(length);
    while (--length >= 0) {
      array[length] = new Vec2();
    }
    return array;
  }
}

class Mat2 {
  constructor(a, b, c, d) {
    this.set(a, b, c, d);
  }
  set(a = 0, b, c, d) {
    if (!b) {
      if (typeof a === 'number') {
        const c = Math.cos(a);
        const s = Math.sin(a);
        this.m00 = c;
        this.m01 = -s;
        this.m10 = s;
        this.m11 = c;
      } else if (a.m00) {
        this.m00 = a.m00;
        this.m01 = a.m01;
        this.m10 = a.m10;
        this.m11 = a.m11;
      }
    } else if (a && b && c && d) {
      this.m00 = a;
      this.m01 = b;
      this.m10 = c;
      this.m11 = d;
    }
    return this;
  }
  absi() {
    return this.abs(this);
  }
  abs(out = new Mat2()) {
    out.m00 = Math.abs(this.m00);
    out.m01 = Math.abs(this.m01);
    out.m10 = Math.abs(this.m10);
    out.m11 = Math.abs(this.m11);
    return this;
  }
  getAxisX(out = new Vec2()) {
    out.x = this.m00;
    out.y = this.m10;
    return out;
  }
  getAxisY(out = new Vec2()) {
    out.x = this.m01;
    out.y = this.m11;
    return out;
  }
  transposei() {
    const t = this.b;
    this.b = this.c;
    this.c = t;
    return this;
  }
  transpose(out = new Mat2()) {
    out.m00 = this.m00;
    out.m01 = this.m10;
    out.m10 = this.m01;
    out.m11 = this.m11;
    return out;
  }
  muli(v) {
    if (typeof v.x === 'number') {
      return this.mul(v.x, v.y, v);
    } else if (v.m00) {
      return this.set(
        this.m00 * v.m00 + this.m01 * v.m10,
        this.m00 * v.m01 + this.m01 * v.m11,
        this.m10 * v.m00 + this.m11 * v.m10,
        this.m10 * v.m01 + this.m11 * v.m11
      );
    }
  }
  mul(x, y, out) {
    // if (typeof x.m00 !== 'number') {
    if (typeof x === 'number') {
      out = out || new Vec2();
      out.x = this.m00 * x + this.m01 * y;
      out.y = this.m10 * x + this.m11 * y;
      return out;
    } else if (typeof x.x === 'number') {
      y = y || new Vec2();
      y.x = this.m00 * x.x + this.m01 * x.y;
      y.y = this.m10 * x.x + this.m11 * x.y;
      return y;
    }
    /* } else {
      out = out || new Mat2();
      out.m00 = this.m00 * x.m00 + this.m01 * x.m10;
      out.m01 = this.m00 * x.m01 + this.m01 * x.m11;
      out.m10 = this.m10 * x.m00 + this.m11 * x.m10;
      out.m11 = this.m10 * x.m01 + this.m11 * x.m11;
      return out;
    } */
  }
}

const ImpulseMath = {
  EPSILON: 0.0001,
  BIAS_RELATIVE: 0.95,
  BIAS_ABSOLUTE: 0.01,
  DT: 1 / 60,
  GRAVITY: new Vec2(0, 50),
  PENETRATION_ALLOWANCE: 0.05,
  PENETRATION_CORRECTION: 0.4,
  clamp(min, max, a) {
    return (a < min ? min : (a > max ? max : a));
  },
  random(min, max) {
    if (String(min).indexOf('.') === -1 && String(max).indexOf('.') === -1) {
      return Math.random() * (max - min) + min;
    } else {
      return Math.floor(Math.random() * (max - min + 1)) + min;
    }
  }
};
ImpulseMath.EPSILON_SQ = ImpulseMath.EPSILON * ImpulseMath.EPSILON;
ImpulseMath.RESTING = ImpulseMath.GRAVITY.mul(ImpulseMath.DT).lengthSq() + ImpulseMath.EPSILON;
ImpulseMath.equal = (a, b) => Math.abs(a - b) <= ImpulseMath.EPSILON;
ImpulseMath.gt = (a, b) => a >= b * ImpulseMath.BIAS_RELATIVE + a * ImpulseMath.BIAS_ABSOLUTE;

class Body {
  constructor(shape, x, y) {
    this.shape = shape;
    this.position = new Vec2(x, y);
    this.velocity = new Vec2(0, 0);
    this.angularVelocity = 0;
    this.torque = 0;
    this.orient = ImpulseMath.random(-Math.PI, Math.PI);
    this.force = new Vec2(0, 0);
    this.staticFriction = 0.5;
    this.dynamicFriction = 0.3;
    this.restitution = 0.2;
    shape.body = this;
    shape.initialize();
  }
  applyForce(f) {
    this.force.addi(f);
  }
  applyImpulse(impulse, contactVector) {
    this.velocity.addsi(impulse, this.invMass);
    this.angularVelocity += this.invInertia * contactVector.cross(impulse);
  }
  setStatic() {
    this.inertia = 0;
    this.invInertia = 0;
    this.mass = 0;
    this.invMass = 0;
  }
  setOrient(radians) {
    this.orient = radians;
    this.shape.setOrient(radians);
  }
}

const Type = Object.freeze({
  Circle: 0,
  Poly: 1,
  Count: 2
});

class Shape {
  constructor() {
    this.u = new Mat2();
  }
  setOrient(radians) {}
}

class Circle extends Shape {
  constructor(r) {
    super();
    this.radius = r;
  }
  clone() {
    return new Circle(this.radius);
  }
  initialize() {
    this.computeMass(1);
  }
  computeMass(density) {
    this.body.mass = Math.PI * this.radius * this.radius * density;
    this.body.invMass = this.body.mass !== 0 ? 1 / this.body.mass : 0;
    this.body.inertia = this.body.mass * this.radius * this.radius;
    this.body.invInertia = this.body.inertia !== 0 ? 1 / this.body.inertia : 0;
    return this;
  }
  getType() {
    return Type.Circle;
  }
}

const MAX_POLY_VERTEX_COUNT = 64;

class Polygon extends Shape {
  constructor(verts, hh) {
    super();
    this.vertices = Vec2.arrayOf(MAX_POLY_VERTEX_COUNT);
    this.normals = Vec2.arrayOf(MAX_POLY_VERTEX_COUNT);
    if (verts instanceof Array) {
      this.set(verts);
    } else if (verts && hh) {
      this.setBox(verts, hh);
    }
  }
  clone() {
    const p = new Polygon();
    p.u.set(this.u);
    for (let i = 0; i < this.vertices.length; ++i) {
      p.vertices[i].set(this.vertices[i]);
      p.normals[i].set(this.normals[i]);
    }
    return p;
  }
  initialize() {
    return this.computeMass(1);
  }
  computeMass(density) {
    const c = new Vec2(0, 0);
    let area = 0;
    let I = 0;
    const k_inv3 = 1 / 3;
    for (let i = 0; i < this.vertices.length; ++i) {
      const p1 = this.vertices[i];
      const p2 = this.vertices[(i + 1) % this.vertices.length];
      const D = p1.crossV(p2);
      const triangleArea = 0.5 * D;
      area += triangleArea;
      const weight = triangleArea * k_inv3;
      c.addsi(p1, weight);
      c.addsi(p2, weight);
      const intx2 = p1.x * p1.x + p2.x * p1.x + p2.x * p2.x;
      const inty2 = p1.y * p1.y + p2.y * p1.y + p2.y * p2.y;
      I += (0.25 * k_inv3 * D) * (intx2 + inty2);
    }
    c.muli(1 / area);
    for (let i = 0; i < this.vertices.length; ++i) {
      this.vertices[i].subi(c);
    }
    this.body.mass = density * area;
    this.body.invMass = this.body.mass !== 0 ? 1 / this.body.mass : 0;
    this.body.inertia = I * density;
    this.body.invInertia = this.body.inertia !== 0 ? 1 / this.body.inertia : 0;
    return this;
  }
  setOrient(radians) {
    this.u.set(radians);
  }
  getType() {
    return Type.Poly;
  }
  setBox(hw, hh) {
    this.vertices = Vec2.arrayOf(4);
    this.vertices[0].set(-hw, -hh);
    this.vertices[1].set(hw, -hh);
    this.vertices[2].set(hw, hh);
    this.vertices[3].set(-hw, hh);
    this.normals = Vec2.arrayOf(4);
    this.normals[0].set(0, -1);
    this.normals[1].set(1, 0);
    this.normals[2].set(0, 1);
    this.normals[3].set(-1, 0);
    return this;
  }
  set(verts) {
    let rightMost = 0;
    let highestXCoord = verts[0].x;
    for (let i = 1; i < verts.length; ++i) {
      const x = verts[i].x;
      if (x > highestXCoord) {
        highestXCoord = x;
        rightMost = i;
      } else if (x === highestXCoord) {
        if (verts[i].y < verts[rightMost].y) {
          rightMost = i;
        }
      }
    }
    const hull = new Array(MAX_POLY_VERTEX_COUNT);
    let outCount = 0;
    let indexHull = rightMost;
    for (;;) {
      hull[outCount] = indexHull;
      let nextHullIndex = 0;
      for (let i = 1; i < verts.length; ++i) {
        if (nextHullIndex === indexHull) {
          nextHullIndex = i;
          continue;
        }
        const e1 = verts[nextHullIndex].clone().subi(verts[hull[outCount]]);
        const e2 = verts[i].clone().subi(verts[hull[outCount]]);
        const c = e1.crossV(e2);
        if (c < 0) {
          nextHullIndex = i;
        }
        if (c === 0 && e2.lengthSq() > e1.lengthSq()) {
          nextHullIndex = i;
        }
      }
      ++outCount;
      indexHull = nextHullIndex;
      if (nextHullIndex === rightMost) {
        this.vertices.length = outCount;
        break;
      }
    }
    for (let i = 0; i < this.vertices.length; ++i) {
      this.vertices[i] = new Vec2(verts[hull[i]]);
    }
    for (let i = 0; i < this.vertices.length; ++i) {
      const face = this.vertices[(i + 1) % this.vertices.length].clone().subi(this.vertices[i]);
      this.normals[i] = new Vec2(face.y, -face.x);
      this.normals[i].normalize();
    }
  }
  getSupport(dir) {
    let bestProjection = -Number.MAX_VALUE;
    let bestVertex = null;
    for (let i = 0; i < this.vertices.length; ++i) {
      const v = this.vertices[i];
      const projection = v.dot(dir);
      if (projection > bestProjection) {
        bestVertex = v;
        bestProjection = projection;
      }
    }
    return bestVertex;
  }
}

const render = (body, ctx) => {
	ctx.strokeStyle = 'white';
  ctx.beginPath();
  if (body.shape.radius) {
    ctx.arc(body.position.x, body.position.y, body.shape.radius, 0, Math.PI * 2, false);
  } else if (body.shape.vertices) {
    for (let i = 0; i < body.shape.vertices.length; ++i) {
      const v = new Vec2(body.shape.vertices[i]);
      body.shape.u.muli(v);
      v.addi(body.position);
      if (i === 0) {
        ctx.moveTo(v.x, v.y);
      } else {
        ctx.lineTo(v.x, v.y);
      }
    }
  }
  ctx.closePath();
  ctx.stroke();
};

const circleCircle = (m, a, b) => {
  const A = a.shape;
  const B = b.shape;
  const normal = b.position.sub(a.position);
  const dist_sqr = normal.lengthSq();
  const radius = A.radius + B.radius;
  if (dist_sqr >= radius * radius) {
    return;
  }
  const distance = Math.sqrt(dist_sqr);
  if (distance === 0) {
    m.penetration = A.radius;
    m.normal.set(1, 0);
    m.contacts[0] = a.position.clone();
  } else {
    m.penetration = radius - distance;
    m.normal.set(normal).divi(distance);
    m.contacts[0] = new Vec2(m.normal).muli(A.radius).addi(a.position);
  }
};

const circlePolygon = (m, a, b) => {
  const A = a.shape;
  const B = b.shape;
  let center = B.u.transpose().muli(a.position.sub(b.position));
  center = center || a.position.sub(b.position);
  let separation = -3.4028235e+38;
  let faceNormal = 0;
  for (let i = 0; i < B.vertices.length; ++i) {
    let s = B.normals[i].dot(center.clone().sub(B.vertices[i]));
    if (s > A.radius) {
      return;
    }
    if (s > separation) {
      separation = s;
      faceNormal = i;
    }
  }
  let v1 = B.vertices[faceNormal];
  let v2 = B.vertices[faceNormal + 1 < B.vertices.length ? faceNormal + 1 : 0];
  if (separation < ImpulseMath.EPSILON) {
    B.u.mul(B.normals[faceNormal], m.normal).negi();
    m.contacts[0] = new Vec2(m.normal).muli(A.radius).addi(a.position);
    m.penetration = A.radius;
    return;
  }
  let dot1 = center.sub(v1).dot(v2.sub(v1));
  let dot2 = center.sub(v2).dot(v1.sub(v2));
  m.penetration = A.radius - separation;
  if (dot1 <= 0) {
    if (center.distanceSq(v1) > A.radius * A.radius) {
      return;
    }
    B.u.muli(m.normal.set(v1).subi(center));
    m.normal.normalize();
    B.u.mul(v1, m.contacts[0] = new Vec2());
    m.contacts[0].addi(b.position);
  } else if (dot2 <= 0) {
    if (center.distanceSq(v2) > A.radius * A.radius) {
      return;
    }
    B.u.muli(m.normal.set(v2).subi(center));
    m.normal.normalize();
    B.u.mul(v2, m.contacts[0] = new Vec2()).addi(b.position);
  } else {
    const n = B.normals[faceNormal];
    if (center.sub(v1).dot(n) > A.radius) {
      return;
    }
    B.u.mul(n, m.normal);
    m.normal.negi();
    const norm = n.clone();
    norm.negi();
    norm.rotate(b.orient);
    m.contacts[0] = new Vec2(a.position).addsi(norm, A.radius);
  }
};

