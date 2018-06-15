var atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'];
var weights = [1.008, 4.002, 6.941, 9.012, 10.814, 12.011, 14.007, 15.999, 18.998, 20.179, 22.989, 24.305, 26.981, 28.085, 30.973, 32.060, 35.450, 39.948, 39.098, 40.078, 44.955, 47.867, 50.941, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.380, 69.723, 72.630, 74.921, 78.971, 79.904, 83.798, 85.468, 87.620, 88.905, 91.224, 92.906, 95.950, 98, 101.070, 102.906, 106.420, 107.868, 112.414, 114.818, 118.710, 121.760, 127.600, 126.90447, 131.293, 132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 145, 150.36, 151.964, 157.25, 158.92535, 162.500, 164.93033, 167.259, 168.93422, 173.045, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966569, 200.592, 204.38, 207.2, 208.98040, 209, 210, 222, 223, 226, 227, 232.0377, 231.03588, 238.02891, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 266, 267, 268, 269, 270, 277, 278, 281, 282, 285, 286, 289, 290, 293, 294, 294];

var regex1 = /[A-Z][^A-Z]*/g;
var last = "H2O+CO2=C6H12O6+O2";

function molecularWeight(molecule) {
  var tWeight = 0;
  if (molecule.match(regex1) === null) {
    return "Invalid";
  }
  molecule.match(regex1).forEach(function(mPart) {
    if (!(/\d+$/.test(mPart))) {
      mPart += "1";
    }
    tWeight += parseFloat(mPart.match(/\d+$/)) * weights[atoms.indexOf(mPart.match(/[A-Z][a-z]*/)[0])];
  });
  return tWeight;
}

function formulaRatio(equation) {
  var molec = [];
  var weight = [];
  var ratio = [];
  var signOrder = [0];
  var totalBS = 0;
  if (equation.match(/[\+\=]/g) === null) {
    return "Invalid";
  }
  equation.match(/[\+\=]/g).forEach(function(sign) {
    if (sign === "+" && totalBS === 0) {
      signOrder.push(0);
    } else {
      signOrder.push(1);
      totalBS = 1;
    }
  });
  equation.match(/[^+=]+/g).forEach(function(molecule) {
    molec.push(molecule);
    if (molecularWeight(molecule.match(/^[\d]*(.*)/)[1]) === "Invalid") {
      return "Invalid 2";
    }
    weight.push(molecularWeight(molecule.match(/^[\d]*(.*)/)[1]));
    if (!(/^[\d]/.test(molecule))) {
      molecule = "1" + molecule;
    }
    ratio.push((parseInt(molecule.match(/^[\d]*/))) * molecularWeight(molecule.match(/^[\d]*(.*)/)[1]));
  });
  return [molec, weight, ratio, equation.match(/[+=]/g).indexOf("="), signOrder];
}

function solveRatio(ratio, known, knownIn) {
  var ratioNew = [];
  var x = known / ratio[knownIn];
  for (var i = 0; i < ratio.length; i++) {
    ratioNew.push(ratio[i] * x);
  }
  return ratioNew;
}

function limitingReagent(matrix, weights, places) {
  var m = (weights[1] * matrix[2][places[0]] / matrix[2][places[1]]);
  var n = (weights[0] * matrix[2][places[1]] / matrix[2][places[0]]);
  console.log(m, n);
  if ((matrix[2][places[0]] / weights[0]) / (matrix[2][places[1]] / weights[1]) > 0.9 && (matrix[2][places[0]] / weights[0]) / (matrix[2][places[1]] / weights[1]) < 1.1) {
    return "balanced to a reasonable amount";
  }
  if (matrix[2][places[0]] / weights[0] > matrix[2][places[1]] / weights[1]) {
    return "reactant " + matrix[0][places[0]] + " is the LLimiting reagent. You need " + m;
  } else {
    return "reactant " + matrix[0][places[1]] + " is the Limiting reagent. You need " + n;
  }
}

function reaKtant(matrix, weights, places) {
  var K = (matrix[2][places[0]] / matrix[2][places[1]]) / (weights[0] / weights[1]);
  if (K > 1) {
    return "Check your numbers and try again!";
  } else {
    return "The reactant efficiancy " + String(parseInt(K * 100000) / 1000) + "%";
  }
}

function Tester(num, ber, aa) {
  console.log(aa);
  aa.backgroundColor = '#f00';
  var matrix = formulaRatio(last);
  //var activeP = 0;
  //var activeR = 0;
  var storage = "";
  var values = [];
  var places = [];
  var known = document.getElementById("cInp" + String(num)).value;
  console.log(solveRatio(matrix[2], known, num));
  var i = 0;
  /*solveRatio(matrix[2], known, num).forEach(function(weight) {
  	var a = document.getElementById("cOup" + String(i));
  	a.innerHTML = parseInt(weight * 1000) / 1000;
  	i++;
  });*/
  for (var j = 0; j < ber; j++) {
    if (document.getElementById("cInp" + String(j)).value !== "") {
      storage += String(matrix[4][j]);
      values.push(document.getElementById("cInp" + String(j)).value);
      places.push(j);
    }
  }
  if (storage === "00") {
    document.getElementById("FOut").innerHTML = limitingReagent(matrix, values, places);
  } else if (storage === "01") {
    document.getElementById("FOut").innerHTML = reaKtant(matrix, values, places);
  } else {
    document.getElementById("FOut").innerHTML = "Please enter in another value";
  }
  console.log(storage, values);
  if (places.length === 1) {
    solveRatio(matrix[2], document.getElementById("cInp" + places[0]).value, places[0]).forEach(function(weight) {
      var a = document.getElementById("cOup" + String(i));
      a.innerHTML = parseInt(weight * 1000) / 1000;
      i++;
    });
  } else {
    for (var k = 0; k < matrix[2].length; k++) {
      document.getElementById("cOup" + String(k)).innerHTML = "NULL";
    }
  }
}

function makeTable(matrix) {
  var thing = document.getElementById("box");
  while (thing.children.length > 1) {
    thing.removeChild(thing.lastChild);
  }
  for (var i = 0; i < matrix[0].length; i++) {
    var box = document.createElement("div");

    var a = document.createElement("p");
    a.innerHTML = '\\(' + matrix[0][i].replace(/\D+/g, '\\text{$&}').replace(/\D(\d+)/g, '}_{$1}') + '\\)';
    a.classList.add("white");
    box.appendChild(a);

    var b = document.createElement("p");
    b.innerHTML = parseInt(matrix[1][i] * 1000) / 1000;
    box.appendChild(b);

    var c = document.createElement("p");
    c.innerHTML = parseInt(matrix[2][i] * 1000) / 1000;
    box.appendChild(c);

    var d = document.createElement("input");
    d.setAttribute("id", "cInp" + String(i));
    d.setAttribute("onKeyUp", "Tester(" + String(i) + "," + String(matrix[0].length) + "," + String("cInp" + String(i)) + ")");
    box.appendChild(d);

    var f = document.createElement("p");
    f.setAttribute("id", "cOup" + String(i));
    f.innerHTML = "NULL";
    box.appendChild(f);

    thing.appendChild(box);

    if (i < matrix[0].length - 1) {
      box = document.createElement("div");

      var e = document.createElement("p");
      e.innerHTML = i === matrix[3] ? "\\(=\\)" : "\\(+\\)";
      e.classList.add("white");
      e.classList.add("sign");
      box.appendChild(e);

      thing.appendChild(box);
    }
  }
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}

makeTable(formulaRatio(last));

function invertMatrix(M) {
  //I didnt make this function. I tried. I failed. SO... you get this code
  //enjoy.
  //maybe I can try to make my own later.
  // source at http://blog.acipo.com/matrix-inversion-in-javascript/
  //just this function
  //i did the rest
  //just to be clear


  //AAAAAAAAAAAAAAAAACCCCCCKKKKK
  //WHAT IS THIS CODE
  //ITS IN SPACES WHATTTTTTT
  //ACACACACACACKCKAKCKAKA
  if (M.length !== M[0].length) {
    return;
  }
  var i = 0,
    ii = 0,
    j = 0,
    dim = M.length,
    e = 0,
    t = 0;
  var I = [],
    C = [];
  for (i = 0; i < dim; i += 1) {
    I[I.length] = [];
    C[C.length] = [];
    for (j = 0; j < dim; j += 1) {
      if (i == j) {
        I[i][j] = 1;
      } else {
        I[i][j] = 0;
      }
      C[i][j] = M[i][j];
    }
  }
  for (i = 0; i < dim; i += 1) {
    e = C[i][i];
    if (e == 0) {
      for (ii = i + 1; ii < dim; ii += 1) {
        if (C[ii][i] != 0) {
          for (j = 0; j < dim; j++) {
            e = C[i][j];
            C[i][j] = C[ii][j];
            C[ii][j] = e;
            e = I[i][j];
            I[i][j] = I[ii][j];
            I[ii][j] = e;
          }
          break;
        }
      }
      e = C[i][i];
      if (e == 0) {
        return
      }
    }
    for (j = 0; j < dim; j++) {
      C[i][j] = C[i][j] / e;
      I[i][j] = I[i][j] / e;
    }
    for (ii = 0; ii < dim; ii++) {
      if (ii == i) {
        continue;
      }
      e = C[ii][i];
      for (j = 0; j < dim; j++) {
        C[ii][j] -= e * C[i][j];
        I[ii][j] -= e * I[i][j];
      }
    }
  }
  return I;
}

function fixReturn(thing) {
  var sum = 0;
  if (thing === null) {
    return 0;
  }
  thing.forEach(function(thing2) {
    //console.log(thing2, thing);
    if (thing2 === "") {
      sum += 1;
    } else {
      sum += parseInt(thing2);
    }
  });
  //console.log(sum);
  return sum;
}

function checkBal(equation) {
  var elems = [];
  var molcs = equation.match(/[A-Z][^+=]*/gm);
  var molcs2 = equation.match(/[^+=]+/gm);
  equation.match(/[A-Z][^+=]*/gm).forEach(function(molecule) {
    console.log(molecule);
    molecule.match(/[A-Z][a-z]?/gm).forEach(function(element) {
      if (elems.indexOf(element) < 0) {
        elems.push(element);
      }
    });
  });
  //console.log(elems, molcs);
  var matrixP = [];
  elems.forEach(function(element) {
    matrixP.push([]);
    molcs.forEach(function(molecule) {
      //console.log("molecule = " + molecule + ", tested element = " + element);
      //console.log("match = " + molecule.match(new RegExp("(?<=" + element + ")[0-9]*","gm")));
      matrixP[matrixP.length - 1].push(fixReturn(molecule.match(new RegExp("(?<=" + element + ")[0-9]*", "gm"))));
    });
  });
  var eqdist = equation.match(/[+=]/gm).indexOf("=") + 1;
  console.log(eqdist);
  var bad = 0
  matrixP.forEach(function(element) {
    var x = element;
    var y = x.splice(eqdist);
    var xSum = 0;
    var ySum = 0;
    var i = 0;
    x.forEach(function(num) {
      xSum += num * fixReturn(molcs2[i].match(/^[\d]*/));
      //console.log(fixReturn(molcs[i].match(/^[\d]*/)));
      i++;
    });
    y.forEach(function(num) {
      ySum += num * fixReturn(molcs2[i].match(/^[\d]*/));
      i++;
      //console.log(fixReturn(molcs[i].match(/^[\d]*/)));
    });
    //console.log(x,y,xSum,ySum);
    console.log(xSum, ySum, x, y);
    console.log(molcs2);
    if (xSum !== ySum) {
      console.log("mrrp");
      bad = 1;
    }
  });
  return bad;
}

function balEq(equation) {
  try {
		if(checkBal(equation)===0){
			return [equation,2];
		}
    var elems = [];
    var molcs = equation.match(/[A-Z][^+=]*/gm);
    equation.match(/[A-Z][^+=]*/gm).forEach(function(molecule) {
      console.log(molecule);
      molecule.match(/[A-Z][a-z]?/gm).forEach(function(element) {
        if (elems.indexOf(element) < 0) {
          elems.push(element);
        }
      });
    });
    console.log(elems, molcs);
    var matrixP = [];
    elems.forEach(function(element) {
      matrixP.push([]);
      molcs.forEach(function(molecule) {
        console.log("molecule = " + molecule + ", tested element = " + element);
        console.log("match = " + molecule.match(new RegExp("(?<=" + element + ")[0-9]*", "gm")));
        matrixP[matrixP.length - 1].push(fixReturn(molecule.match(new RegExp("(?<=" + element + ")[0-9]*", "gm"))));
      });
    });
    console.log(matrixP);
    while (matrixP[0].length > matrixP.length) {
      matrixP.push([]);
      for (var i = 0; i < matrixP[0].length; i++) {
        if (i === matrixP.length - 1) {
          matrixP[matrixP.length - 1].push(1);
        } else {
          matrixP[matrixP.length - 1].push(0);
        }
      }
    }
    var ratio = [];
    var matrix = invertMatrix(matrixP);
    for (var k = 0; k < matrix.length; k++) {
      ratio.push(matrix[k][matrix.length - 1]);
    }
    console.log(ratio);
    var quitValue = 3000;
    var multiplier = 0;
    while (quitValue !== 0) {
      multiplier++;
      quitValue++;
      var bad = 0;
      ratio.forEach(function(numb) {
        if (!(Math.round(numb * 1000 * multiplier) % 1000 === 0)) {
          bad = 1;
        }
      });
      if (bad === 0) {
        break;
      }
    }
    console.log(multiplier);
    var answer = "";
    for (var l = 0; l < ratio.length; l++) {
      ratio[l] = ratio[l] * multiplier;
      if (Math.abs(ratio[l]) !== 1) {
        answer += Math.abs(ratio[l]);
      }
      answer += molcs[l];
      if (ratio.length - 1 > l) {
        if (ratio[l + 1] * ratio[l] > 0) {
          answer += "+";
        } else {
          answer += "=";
        }
      }
    }
    console.log(answer);
		if(checkBal(answer) === 1){
			return ["your+mom",3]
		}
    return [answer, 0];
  } catch (TypeError) {
    return [equation, 1];
  }
}

var init = true;

function thing() {
  var value = document.getElementById("in").value;
  console.log(value);
  if (value !== last || init) {
    init = false;
    last = value;
    var aff = balEq(document.getElementById("in").value);
    last = aff[0];
    var result = formulaRatio(aff[0]);
    if (result === "Invalid") {
      document.getElementById("warning").innerHTML = "Please complete the equation.";
      return;
    } else if (result === "Invalid 2" || !result) {
      document.getElementById("warning").innerHTML = "Invalid input.";
      return;
    }
    makeTable(result);
    if (aff[1] === 0) {
      document.getElementById("warning").innerHTML = "Balancing Successful";
      document.getElementById("warning").style.color = "black";
    }
		if(aff[1] === 1) {
      document.getElementById("warning").innerHTML = "Equation could not be balanced.";
      document.getElementById("warning").style.color = "red";
    }
		if(aff[1] === 2){
			document.getElementById("warning").innerHTML = "Wow. You balanced it.";
      document.getElementById("warning").style.color = "blue";
		}
		if(aff[1] === 3){
			document.getElementById("warning").innerHTML = "Is everything in the correct order?";
      document.getElementById("warning").style.color = "red";
		}
  }
}
thing();
