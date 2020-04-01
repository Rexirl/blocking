self.importScripts('perlin.js', '../../js/three.js', 'Blocking.js');

const setUp = (data) => {
    self.workerIndex = data.workerIndex;
    self.chunkResolution = new THREE.Vector3();
    for (let k of ['x', 'y', 'z']) {
        self.chunkResolution[k] = data.chunkResolution[k];
    }
    self.padding = 1; // padding on one side (so two times for full padding)
    self.padRes = chunkResolution.clone().addScalar(padding * 2);
    self.numBlocks = self.padRes.x * self.padRes.y * self.padRes.z;
    
    // Set up persistent work buffers (ArrayBuffer) for the maximum
    // number of faces possible, to avoid counting in a separate pass below 
    // (rather copy from the work buffer when done)
    self.maxNumFaces = 12 * self.numBlocks; // - actually lower (or perhaps higher, if detail geometry is introduced)
    self.temp = {
        position: {
            components: 3,
            bytesPerComponent: 1,
            type: Uint8Array,
            normalize: false
        },
        normal: {
            components: 3,
            bytesPerComponent: 1,
            type: Int8Array,
            normalize: false
        },
        occlusion: {
            components: 1,
            bytesPerComponent: 1,
            type: Uint8Array,
            normalize: true
        }
    }
    for (let ti in self.temp) {
        const t = self.temp[ti];
        t.array = new ArrayBuffer(1 * self.maxNumFaces * t.components * t.bytesPerComponent);
        t.typedArray = new t.type(t.array);
    }
}

const generateChunkGeometry = (chunk) => {
    const baseNormals = [
        new THREE.Vector3(1, 0, 0),
        new THREE.Vector3(0, 1, 0),
        new THREE.Vector3(0, 0, 1)
    ];
    
    const isEmpty = bi => {
        const type = chunk.blocks[bi];
        return type === BLOCK_TYPES.EMPTY;
    };
    const inPadding = c => {
        for (let k of ['x', 'y', 'z']) {
            if (c[k] < 0 || c[k] >= chunkResolution[k]) return true;
        }
        return false;
    }
    
    chunk.allEmpty = true;
    const t1 = Date.now();
    const position  = self.temp.position.typedArray;
    const normal    = self.temp.normal.typedArray;
    const occlusion = self.temp.occlusion.typedArray;
    
    let vbase = [];
    let basei = 0;
    for (let z = 0; z < 2; ++z)
    for (let y = 0; y < 2; ++y)
    for (let x = 0; x < 2; ++x) {
        vbase.push(new THREE.Vector3(x, y, z));
    }
    
    const dimAdd = {x: 1, y: padRes.x, z: padRes.x * padRes.y};
    let ii = 0;
    let vi = 0, ni = 0, oi = 0;
    chunk.blockIteration((c, bi) => {
        if (inPadding(c)) return;
        const type = chunk.blocks[bi];
        if (isEmpty(bi)) return;
        chunk.allEmpty = false;
        
        // Make faces.
        const faces = 
        [[0, 4, 2], [2, 4, 6], [1, 3, 5], [5, 3, 7],    // left, right
         [0, 1, 4], [4, 1, 5], [3, 2, 6], [3, 6, 7],    // bottom, top
         [1, 0, 3], [3, 0, 2], [4, 5, 7], [4, 7, 6]];   // front, back (no, actually reversed order)
        const makeFace = (i, nk, ns) => {
            i *= 2;
            const n = new THREE.Vector3();
            n[nk] = ns;
            
            const faces = [
                [0, 4, 6, 2], [1, 3, 7, 5], // left, right
                [0, 1, 5, 4], [3, 2, 6, 7], // bottom, top
                [1, 0, 2, 3], [4, 5, 7, 6], // back, front
            ];
            const fv = [];
            fv.length = 4;
            for (let fvi = 0; fvi < fv.length; ++fvi) {
                const v = vbase[faces[i / 2][fvi]].clone();
                
                // Compute occlusion factor:
                let occ = 1.0; 
                const nAdd = ns * dimAdd[nk];
                let totalAdd = nAdd;
                let sides = [];
                for (let k of ['x', 'y', 'z']) {
                    if (k == nk) continue;
                    const addAdd = (v[k] * 2 - 1) * dimAdd[k];
                    sides.push(isEmpty(bi + nAdd + addAdd) ? 0 : 1);
                    totalAdd += addAdd;
                }
                const corner = isEmpty(bi + totalAdd) ? 0 : 1;
                occ = (sides[0] && sides[1]) ? 0 : 3 - (sides[0] + sides[1] + corner);
                occ /= 3;
                
                v.add(c);
                fv[fvi] = {p: v, n: n, o: occ};
            }
            
            // Flip the face split depending consistently on occlusion:
            const flip = fv[0].o + fv[2].o > fv[1].o + fv[3].o;
            for (let ii = 0; ii < 6; ++ii) {
                const ti = (flip ? [0, 1, 2, 2, 3, 0] : [0, 1, 3, 3, 1, 2])[ii];
                const v = fv[ti];
                occlusion[oi++] = v.o * 255;
                for (let k of ['x', 'y', 'z']) {
                    position[vi++] = v.p[k];
                    normal[ni++] = n[k];
                }
            }
        };
        let add = 1;
        for (let d = 0; d < 3; ++d) {
            const k = (['x', 'y', 'z'])[d];
            let n = baseNormals[d];
            if (isEmpty(bi - add)) {
                makeFace(d * 2, k, -1); // lower
            }
            if (isEmpty(bi + add)) {
                makeFace(d * 2 + 1, k, 1); // higher
            }
            add *= padRes[k];
        }
    });
    
    const numFaces = vi / 3;
    if (!numFaces) return;
    
    chunk.mesh = {};
    for (let ti in self.temp) {
        const t = self.temp[ti];
        chunk.mesh[ti] = {components: t.components, normalize: t.normalize};
    }
    
    // Create new arrays of correct sizes and copy to them
    const copyAttributeArray = (t) => {
        return t.array.slice(0, numFaces * t.components * t.bytesPerComponent);
    }
    chunk.positionArray   = copyAttributeArray(temp.position);
    chunk.normalArray     = copyAttributeArray(temp.normal);
    chunk.occlusionArray  = copyAttributeArray(temp.occlusion);
    
    //console.log('Face creation in: ' + (Date.now() - t1));
}

self.onmessage = e => {
    const t1 = Date.now();
    const data = e.data[0];
    if (data.setUp) {
        setUp(data);
        return;
    }
    
    const chunkPos = new THREE.Vector3().copy(data.position);
    const chunkResolution = self.chunkResolution;//new THREE.Vector3().copy(data.chunkResolution);
    
    const padding = self.padding;
    const padRes = self.padRes;
    
    const numBlocks = self.numBlocks;
    let blocksArray = new ArrayBuffer(2 * numBlocks); // will use type Uint16
    
    const chunk = {};
    // f: function with parameters ({x: ix, y: iy, z: iz}, bi)
    chunk.blockIteration = (f) => {
        let c = new THREE.Vector3();
        let bi = 0;
        for (c.z = -padding; c.z < chunkResolution.z + padding; ++c.z)
        for (c.y = -padding; c.y < chunkResolution.y + padding; ++c.y)
        for (c.x = -padding; c.x < chunkResolution.x + padding; ++c.x) {
            f(c, bi);
            ++bi;
        }
    }
    
    const generateChunk = (chunk) => {
        const id = chunkPos.clone().divide(chunkResolution);
        chunk.id = id;
        chunk.position = chunkPos;
        chunk.key = JSON.stringify(id);
        chunk.blocks = new Uint16Array(blocksArray);
        
        chunk.blockIteration((c, bi) => {
            const ac = chunk.position.clone().add(c);
            let density = 0.0;
            
            { // Shapes (to be improved):
                
                const center = new THREE.Vector3(9990, 0, 0);
                const radius = 10000;
                const normal = new THREE.Vector3(0, 1, 0);
                const proj = ac.clone().sub(center).sub(normal.clone().multiplyScalar(ac.dot(normal)));
                
                let noiseValue = 0.0;
                { // Noise:
                    let v = 0.0;
                    let c = ac.clone();
                    let a = 0.5;
                    const numOctaves = 10; // - doubles block generation time at 10 compared to 1. To do: custom noise that evaluates fewer times at lower detail, if possible
                    for (let i = 0; i < numOctaves; ++i) {
                        v += a * noise.simplex3(c.x, c.y, c.z);
                        c.divideScalar(2.0);
                        a *= 2;
                    }
                    //if (v > 0.0) density = 1.0;
                    //dist += v;
                    noiseValue = v;
                }
                { // Ground:
                    //if (ac.y < 2) density = 1.0;
                    let dist = Math.max(ac.y, -2) - 2;
                    dist += noiseValue * 0.1;
                    
                    
                    if (dist < 0.0) density = 1.0;
                }
                {// "Stairs"/"Amphitheatre":
                    const stairs = () => {
                        const height = 80;
                        const dist = ac.distanceTo(center);
                        if (dist < radius) return false;
                        
                        if (ac.y > dist - radius) return false;
                        //if (ac.x > -10) return false;
                        //if (-10 - ac.x - ac.y <= 0) return false;
                        
                        if (ac.y > height + noiseValue) return false;
                        if (ac.y < 0.0 + noiseValue) return false;
                        return true;
                    };
                    if (stairs()) density = 1.0;
                }
                
                const getTowerBase = () => {
                    const towerPosRadius = 9990;
                    const numTowers = 729.0;
                    
                    const np = proj.clone().normalize();
                    let ti = Math.round(Math.atan2(np.z, np.x) / Math.PI * numTowers);
                    const angle = ti / numTowers * Math.PI;
                    let base = new THREE.Vector3(Math.cos(angle), 0, Math.sin(angle)).multiplyScalar(towerPosRadius).add(center);
                    base.x = Math.floor(base.x);
                    base.z = Math.floor(base.z);
                    return {ti: ti, base: base};
                }
                const towerBase = getTowerBase();
                
                {// Subterranean entrance(s):
                    const ti = towerBase.ti;
                    //if (ti == 729 || ti == 728) density = 0;
                    const cutOut = () => {
                        if (proj.z <= 5 || proj.z >= 38) return false;
                        if (proj.length() < radius + 6) return false;
                        //if (proj.length() > radius + 40) return false;
                        if (proj.x + radius - ac.y + 128 < 0) return false;
                        if (proj.x + radius - ac.y + 12 > 0) return false;
                        //if (Math.abs(proj.x) - radius - ac.y < 0) return false;
                        return true;
                    }
                    if (ac.y > -10000 && cutOut()) {
                        density = 0;
                        if (noiseValue - Math.max(-20, ac.y) > 0 && ac.y < 0) density = 1;
                    }
                }
                {// Towers:
                    const tower = () => {
                        const base = towerBase.base;
                        const height = 64;
                        let width = 7;
                        let widthZ = 4;
                        if (Math.abs(ac.y % 8) == 7) {
                            width -= 1;
                            widthZ -= 1;
                        }
                        
                        if (ac.z <= base.z - widthZ) return false;
                        if (ac.z >= base.z + widthZ) return false;
                        if (ac.x <= base.x - width ) return false;
                        if (ac.x >= base.x + width ) return false;
                        if (ac.y >= base.y + height) return false;
                        return true;
                    }
                    if (tower()) density = 1.0;
                    //if (ac.z > -20 && ac.z < -10 && ac.x > -5 && ac.x < 5 && ac.y < 64) density = 1.0;
                }
            }
            
            let blockType = density > 0 ? BLOCK_TYPES.FILLED : BLOCK_TYPES.EMPTY;
            chunk.blocks[bi] = blockType;
        });
    }
    generateChunk(chunk);
    const t2 = Date.now();
    
    //let positionArray, normalArray, occlusionArray; // attributes
    
    generateChunkGeometry(chunk);
    let positionArray = chunk.positionArray,
        normalArray = chunk.normalArray;
        occlusionArray = chunk.occlusionArray;
    let transferables = [positionArray, normalArray, occlusionArray].filter(e => e != undefined);
    for (let k of ['blocks', 'blockIteration', 'positionArray', 'normalArray', 'occlusionArray']) {
        chunk[k] = undefined; // for cloning/transfer
    }
    //chunk.blocks = undefined; // for cloning/transfer
    
    if (chunk.allEmpty) blocksArray = undefined;
    else transferables.push(blocksArray);
    
    //console.log('Generation in ' + (Date.now() - t1) + ' (blocks in ' + (t2 - t1) + ')');
    
    postMessage({
            workerIndex: self.workerIndex,
            chunk: chunk, 
            blocks: blocksArray,
            position: positionArray,
            normal: normalArray,
            occlusion: occlusionArray
        },
        transferables);
};
