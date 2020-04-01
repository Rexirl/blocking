/* Blocking code */

const BLOCK_TYPES = {
    INVALID : 0,
    EMPTY   : 1,
    // - to do: expand/change
    FILLED  : 2
};

class BlockWorld {
    constructor(scene, mat) {
        this.frame           = 0;
        this.lastChunkPosition = {};
        
        this.posToChunk      = {};
        //this.chunks          = [];
        this.scene           = scene;
        this.viewDistance    = 1; // in chunks (to do: try to gain performance for 8, at least)
        this.chunkResolution = new THREE.Vector3(32, 64, 32);   // in blocks
        
        this.materials = {};
        this.materials.wireframe = new THREE.MeshBasicMaterial({
            color: 0xffffff,
            wireframe: true
        });
        this.materials.basic = new THREE.MeshLambertMaterial({ color: 0xF7F7F7 });
        this.material = this.materials.basic;
        //this.material = this.materials.wireframe;
        this.material = mat;
        
        this.totalGenerationTime = 0.0;
        this.numGenerated = this.numNoMesh = this.numRemoved = 0;
        
        this.generators = 4; // number of workers we're allowed to use
        this.workers = [];
        // - 8 workers do make a significant difference (at least on an AMD 2700X)
        this.workers.length = 8; // - somewhat arbitrary
        for (let i = 0; i < this.workers.length; ++i) {
            const worker = new Worker('js/Generator.js');
            this.workers[i] = {worker: worker, num: 0};
            
            worker.postMessage([{
                setUp: true,
                workerIndex: i,
                chunkResolution: this.chunkResolution,
            }]);
            
            worker.onmessage = e => {
                ++this.numGenerated;
                const i = e.data.workerIndex;
                --this.workers[i].num;
                const chunk = e.data.chunk;
                chunk.blocks = e.data.blocks ? new Uint16Array(e.data.blocks) : undefined;
                chunk.loaded = true;
                this.posToChunk[chunk.key] = chunk;
                
                if (e.data.position && e.data.position.byteLength) {
                    const g = new THREE.BufferGeometry();
                    // - to do: re-automate this, despite annoying ArrayBuffer requirements for transfer
                    /*for (let name in chunk.mesh) {
                        const a = chunk.mesh[name];
                        g.addAttribute(name, new THREE.BufferAttribute(e.data[name], a.components, !!a.normalize));
                    }*/
                    g.addAttribute('position', new THREE.BufferAttribute(new Uint8Array(e.data.position), 3));
                    g.addAttribute('normal', new THREE.BufferAttribute(new Int8Array(e.data.normal), 3));
                    g.addAttribute('occlusion', new THREE.BufferAttribute(new Uint8Array(e.data.occlusion), 1, true));
                    chunk.mesh = new THREE.Mesh(g, this.material);
                    this.scene.add(chunk.mesh);
                } else {
                    chunk.mesh = undefined;
                    ++this.numNoMesh;
                }
            };
        }
        
        this.maxPerWorker = 4; // - to tune
        this.ageLimit = 40;//20; // frames before dropping disused chunks (must be higher than view distance because of "shell updating" in partitions
    }
    
    setWireframe(w) {
        /*this.material = w ? this.materials.wireframe : this.materials.basic;
        for (const chunk of this.chunks) {
            chunk.mesh.material = this.material;
            chunk.mesh.material.wireframe = w;
        }*/
    }
    
    // - to do: write correct ray intersect (with always bounded distance)
    // Returns distance to closest non-empty block in d direction, or zero if inside non-empty already
    /*distance(p, d) {
        let minDistance = Infinity;
        for (const k in this.posToChunk) {
            const chunk = this.posToChunk[k];
            if (!chunk.loaded) continue;
            const min = new THREE.Vector3(chunk.position.x, chunk.position.y, chunk.position.z);
            const max = min.clone().add(this.chunkResolution);
            const box = new THREE.Box3(min, max);
            // - to do: don't just check on containment, and give proper distance
            if (box.containsPoint(p)) {
                const rel = p.clone().sub(min).divide(max.sub(min));
                rel.multiply(this.chunkResolution);
                const bi = Math.floor(rel.x) + this.chunkResolution.x * (Math.floor(rel.y) + this.chunkResolution.y * Math.floor(rel.z));
                if (!chunk.blocks[bi]) continue;
                
                const type = chunk.blocks[bi];
                if (type === BLOCK_TYPES.INVALID || type === BLOCK_TYPES.EMPTY) continue;
                return 0.0;
            }
        }
        return minDistance;
    }*/
    
    // p0: THREE.Vector3 ray origin
    // p1: THREE.Vector3 ray end
    // - to do: allow checking a given height (maybe even more general shapes later on)
    // (or maybe that's better dealt with by intersecting a couple of rays)
    // Returns only the first hit (maybe allow for configurability there).
    // Ignores first hit if the origin is inside the block.
    rayIntersect(p0, p1) {
        // - to do, one day (or night): take care to handle numerical imprecision
        
        const result = {distance: Infinity, block: undefined};
        if (p0.equals(p1)) return result; // no hit
        const positionToBlock = p => p.clone().floor();
        const block0 = positionToBlock(p0), block1 = positionToBlock(p1);
        const blockToChunk = i => i.clone().divide(this.chunkResolution).floor();
        const chunk0 = blockToChunk(block0.clone().min(block1)), 
              chunk1 = blockToChunk(block1.clone().max(block1));
        
        const zero3 = new THREE.Vector3();
        const resSub1 = this.chunkResolution.clone().subScalar(1);
        const padding = 1; // - to do: make this a global/static constant somewhere
        const padRes = this.chunkResolution.clone().addScalar(padding * 2);
        
        const d = p1.clone().sub(p0);
        const maxDistance = d.length();
        d.normalize();
        const dinv = new THREE.Vector3(1, 1, 1).divide(d);
        
        // - to do: look into optimising this (3D line drawing through the chunks)
        // (i.e. traverse from current rather than looping through bounds)
        
        let numChunks = 0, numChunksHit = 0, numBlocks = 0, numBlocksHit = 0;
        
        for (let cz = chunk0.z; cz <= chunk1.z; ++cz)
        for (let cy = chunk0.y; cy <= chunk1.y; ++cy)
        for (let cx = chunk0.x; cx <= chunk1.x; ++cx) {
            ++numChunks;
            const chunkId = new THREE.Vector3(cx, cy, cz);
            const chunkPos = chunkId.clone().multiply(this.chunkResolution);
            // - to do: skip chunk if not intersected by ray
            const chunk = this.chunkFromId(chunkId);
            if (chunk && !chunk.blocks) continue; // empty chunk
            
            const aabb = this.computeChunkBox(chunkId);
            const aabbVsRay = box => {
                let t1 = (box.min.x - p0.x) * dinv.x;
                let t2 = (box.max.x - p0.x) * dinv.x;
                let tmin = Math.min(t1, t2);
                let tmax = Math.max(t1, t2);
                for (let k of ['y', 'z']) { // x is handled already
                    t1 = (box.min[k] - p0[k]) * dinv[k];
                    t2 = (box.max[k] - p0[k]) * dinv[k];
                    tmin = Math.max(tmin, Math.min(t1, t2));
                    tmax = Math.min(tmax, Math.max(t1, t2));
                }
                return {
                    tmin: tmin, tmax: tmax,
                    hit: tmax >= tmin//tmax >= Math.max(tmin, 0)
                };
            }
            const i = aabbVsRay(aabb);
            if (!i.hit) continue; // no intersection
            ++numChunksHit;
            
            // - this check isn't wholly accurate (when chunks in the search space aren't loaded
            // but there's space in loaded chunks that's closer than any in the unloaded could be)
            if (!chunk || !chunk.loaded) {
                return undefined; // - maybe don't do this but rather skip the chunk
            }
            
            
            // - this can be optimised greatly (essentially 3D line-drawing)
            // but let's start off with something simple that works for short distances
            const hit0 = d.clone().multiplyScalar(i.tmin).add(p0);
            const hit1 = d.clone().multiplyScalar(i.tmax).add(p0);
            
            const local0 = positionToBlock(hit0).sub(chunkPos).clamp(zero3, resSub1);
            const local1 = positionToBlock(hit1).sub(chunkPos).clamp(zero3, resSub1);
            const chunkBlock0 = local0.clone().min(local1);
            const chunkBlock1 = local0.clone().max(local1);
            
            // - to do: search closest first for early exit
            for (let bz = chunkBlock0.z; bz <= chunkBlock1.z; ++bz)
            for (let by = chunkBlock0.y; by <= chunkBlock1.y; ++by)
            for (let bx = chunkBlock0.x; bx <= chunkBlock1.x; ++bx) {
                ++numBlocks;
                const bi = (bx + padding) + padRes.x * ((by + padding) + padRes.y * (bz + padding));
                const isEmpty = bi => {
                    const type = chunk.blocks[bi];
                    return type === BLOCK_TYPES.INVALID || type === BLOCK_TYPES.EMPTY;
                };
                if (isEmpty(bi)) continue;
                const blockId = new THREE.Vector3(bx, by, bz).add(chunkPos);
                const aabb = new THREE.Box3(blockId, new THREE.Vector3(1, 1, 1).add(blockId));
                const i = aabbVsRay(aabb);
                if (!i.hit || i.tmin < 0) continue;
                ++numBlocksHit;
                if (i.tmin < result.distance) {
                    result.distance = i.tmin;
                    result.block = blockId;
                }
            }
        }
        if (result.distance <= maxDistance) {
            result.position = d.clone().multiplyScalar(result.distance).add(p0);
            return result;
        }
        return {distance: Infinity, block: undefined};
    }
    
    // - to do: generalise (distance checked)
    // (maybe use sphere vs AABB tests - specify a radius)
    // Get vector to closest empty block or undefined if the closest empty block
    // is more than 4 blocks away.
    // - to do: special value for being inside an unloaded chunk)
    // - to do: consider object/player height
    findNear(p) {
        // - maybe to try: precalculating spaces, compactly, somehow
        // (and for fast checking performance)
        
        // At most 8 chunks have to be checked (so long as max distance doesn't exceed the chunk size).
        const pBlock = p.clone().floor();
        const distance = 4; // - to do: make configurable, probably
        const blockMin = pBlock.clone().subScalar(distance);
        const blockMax = pBlock.clone().addScalar(distance);
        
        const chunkMin = blockMin.clone().divide(this.chunkResolution).floor();
        const chunkMax = blockMax.clone().divide(this.chunkResolution).floor();
        const zero3 = new THREE.Vector3();
        const resSub1 = this.chunkResolution.clone().subScalar(1);
        
        const padding = 1; // - to do: make this a global/static constant somewhere
        const padRes = this.chunkResolution.clone().addScalar(padding * 2);
        
        let vector = new THREE.Vector3();
        const result = {
            empty: { distance: Infinity },
            nonEmpty: { distance: Infinity }
        };
        
        for (let cz = chunkMin.z; cz <= chunkMax.z; ++cz)
        for (let cy = chunkMin.y; cy <= chunkMax.y; ++cy)
        for (let cx = chunkMin.x; cx <= chunkMax.x; ++cx) {
            const chunkId = new THREE.Vector3(cx, cy, cz);
            const chunkPos = chunkId.clone().multiply(this.chunkResolution);
            const chunk = this.chunkFromId(chunkId);
            // - this check isn't wholly accurate (when chunks in the search space aren't loaded
            // but there's space in loaded chunks that's closer than any in the unloaded could be)
            if (!chunk || !chunk.loaded) {
                return undefined;
            }
            
            // - to do: handle wholly empty chunks without looping inside
            const chunkEmpty = !chunk.mesh;
            
            const chunkBlockMin = blockMin.clone().sub(chunkPos).clamp(zero3, resSub1);
            const chunkBlockMax = blockMax.clone().sub(chunkPos).clamp(zero3, resSub1);
            
            // - to do: search closest first for early exit
            for (let bz = chunkBlockMin.z; bz <= chunkBlockMax.z; ++bz)
            for (let by = chunkBlockMin.y; by <= chunkBlockMax.y; ++by)
            for (let bx = chunkBlockMin.x; bx <= chunkBlockMax.x; ++bx) {
                // - to do: check for a space of at least some given height (e.g. 2 blocks for players)
                
                const bi = (bx + padding) + padRes.x * ((by + padding) + padRes.y * (bz + padding));
                const isEmpty = bi => {
                    const type = chunk.blocks[bi];
                    return type === BLOCK_TYPES.INVALID || type === BLOCK_TYPES.EMPTY;
                };
                const offset = 0.5;
                const blockPosI = new THREE.Vector3(bx, by, bz).add(chunkPos);
                const blockPos = blockPosI.clone().addScalar(offset);
                const dest = chunkEmpty || isEmpty(bi) ? result.empty : result.nonEmpty;
                const dist = blockPos.distanceTo(p);
                if (dist < dest.distance) {
                    dest.distance = dist;
                    dest.vector = blockPos.sub(p);
                    dest.blockId = blockPosI;
                    dest.blockOffset = blockPosI.clone().sub(pBlock);
                }
            }
        }
        
        return result;
        
        /*const chunkPos = this.positionToChunkId(p);
        const chunkWorldPos = chunkPos.clone().multiply(this.chunkResolution);
        const chunkAdd = new THREE.Vector3();
        for (let k of ['x', 'y', 'z']) {
            chunkAdd[k] = Math.sign(p[k] - chunkWorldPos[k] - this.chunkResolution[k] / 2);
        }
        
        for (let k of ['x', 'y', 'z'])
        for (let i = 0; i < 2; ++i) {
            const cp = chunkPos.clone();
            cp[k] += i * chunkAdd[k];
            const chunk = this.chunkFromPosition(cp);
            if (!chunk || !chunk.loaded) return undefined; // can't determine result
            // - to do: check blocks
        }*/
        // - to do: return Infinity if no emptiness was found, or correct vector otherwise
    }
    
    positionToChunkId(p) {
        const cp = p.clone().divide(this.chunkResolution);
        return new THREE.Vector3(Math.floor(cp.x), Math.floor(cp.y), Math.floor(cp.z));
    }
    
    chunkFromId(chunkId) {
        return this.posToChunk[JSON.stringify(chunkId)];
    }
    
    update(camera, p) {
        this.generate(camera, p);
        const position = new THREE.Vector3(); // - testing
        //const position = p;
        
        // Recompute and set relative chunk positions
        // (to avoid floating-point imprecision in far-away chunks)
        let t1 = Date.now();
        for (const chunkPos in this.posToChunk) {
            const chunk = this.posToChunk[chunkPos];
            if (!chunk.mesh) continue;
            chunk.mesh.position.setX(chunk.position.x - position.x);
            chunk.mesh.position.setY(chunk.position.y - position.y);
            chunk.mesh.position.setZ(chunk.position.z - position.z);
        }
        
        /*{
            if (!this.dontShowTooOften) this.dontShowTooOften = 1;
            if (this.dontShowTooOften % 60 == 0)
            console.log(`Position update in ${Date.now() - t1}`);
        }*/
    }
    
    // chunkId: position in chunk coordinates (not yet multiplied by chunk block size)
    computeChunkBox(chunkId) {
        const min = chunkId.clone().multiply(this.chunkResolution);
        return new THREE.Box3(min, min.clone().add(this.chunkResolution));
    }
    
    setUpUpdateOrder() {
        this.viewDistance  = Math.floor(this.viewDistance);
        this.viewDistanceY = Math.floor(this.viewDistance / (this.chunkResolution.y / this.chunkResolution.x));
        if (this.viewDistanceY < 1) this.viewDistanceY = 1;
        if (this.viewDistance === this.lastViewDistance) return;
        this.lastViewDistance = this.viewDistance;
        
        let numDistant = 0;
        
        // List of relative-from-camera chunk coordinates to update for successive view distances:
        const maxDist2 = Math.pow(this.chunkResolution.x * (this.viewDistance + 1), 2);
        const range = this.viewDistance;
        const rangeY = this.viewDistanceY;
        const chunkOrder = this.chunkOrder = [];
        for (let iz = -range;  iz <= range; ++iz)
        for (let iy = -rangeY; iy <= rangeY; ++iy)
        for (let ix = -range;  ix <= range; ++ix) {
            const chunkPos = new THREE.Vector3(ix, iy, iz);
            const box = this.computeChunkBox(chunkPos);
            const center = box.max.clone().add(box.min).divideScalar(2.0);
            // - to do: solve inaccuracies (compute more precisely)
            if (center.clone().lengthSq() > maxDist2) { // - inaccurate for y (low view distance especially)
                ++numDistant;
                continue; // - too distant
            }
            chunkOrder.push({position: chunkPos, distance: center.length()});
        }
        chunkOrder.sort((a, b) => a.distance - b.distance);
    }
    
    generate(camera, position) {
        /*{
            const closest = this.vectorToEmpty(position);
            if (closest && this.dbgVectorToEmpty % 60 === 0)
            console.log('Closest distance: ' + closest.distance + ', vector: (', 
                + closest.vector.x + ', ' + closest.vector.y + ', ' + closest.vector.z + ')');
        }*/
        this.setUpUpdateOrder();
        const t1 = Date.now();
        
        let frustum = new THREE.Frustum();
        frustum.setFromMatrix(new THREE.Matrix4().multiply(camera.projectionMatrix).multiply(camera.matrixWorldInverse));
        
        let posChunk = new THREE.Vector3();
        let oldEquality = true;
        for (let k of ['x', 'y', 'z']) {
            posChunk[k] = Math.floor(position[k] / this.chunkResolution[k]);
            oldEquality = posChunk[k] === this.lastChunkPosition[k] && oldEquality;
        }
        //if (oldEquality) return; // no need for an update (this'll change with frustum consideration)
        
        this.lastChunkPosition = posChunk;
        
        let numVisited = 0, numDistant = 0;
        const maxDist2 = Math.pow(this.chunkResolution.x * (this.viewDistance + 1), 2);
        
        const range = this.viewDistance;
        const rangeY = this.viewDistanceY;
        
        const chunksExtent = new THREE.Vector3(range, rangeY, range);
        const minChunkPos = posChunk.clone().sub(chunksExtent);
        const maxChunkPos = posChunk.clone().add(chunksExtent);
        

        const ipart = 1024;
        const inum = Math.ceil(this.chunkOrder.length / ipart);
        let i0 = ipart * (this.frame % inum);
        for (let ci = i0; ci < Math.min(this.chunkOrder.length, i0 + ipart); ++ci) {
            ++numVisited;
            const chunkPos = this.chunkOrder[ci].position.clone().add(posChunk);
            const box = this.computeChunkBox(chunkPos);
            const center = box.max.clone().add(box.min).divideScalar(2.0);
            
            // - to do: perhaps special case for a few chunks right next to the camera?
            //if (!frustum.intersectsBox(box)) continue;
            
            const chunkKey = JSON.stringify(chunkPos); // - can we get rid of this? Probably (multi-dimensional array instead?)
            if (this.posToChunk[chunkKey]) {
                this.posToChunk[chunkKey].frame = this.frame;
                continue;
            } else if (ci < 27) { // - quick and untuned way of always loading the closest chunks
            } else if (!frustum.intersectsBox(box)) {
                continue;
            }
            
            // Asynchronous loading with web workers:
            const findWorker = () => {
                let leastNum = Infinity;
                let worker;
                for (let i = 0; i < Math.min(this.workers.length, this.generators); ++i) {
                    let w = this.workers[i];
                    if (w.num < leastNum) {
                        leastNum = w.num;
                        worker = w;
                    }
                }
                return worker;
            }
            const worker = findWorker();
            if (worker.num >= this.maxPerWorker) continue;
            ++worker.num;
            this.posToChunk[chunkKey] = {};
            chunkPos.multiply(this.chunkResolution);
            let workerIndex = 0
            for (; workerIndex < this.workers.length; ++workerIndex) {
                if (this.workers[workerIndex] === worker) break;
            }
            worker.worker.postMessage([{
                position: chunkPos,
            }]);
        }
        const t2 = Date.now();
        
        // Remove unused chunks:
        for (const chunkPos in this.posToChunk) {
            const chunk = this.posToChunk[chunkPos];
            if (!chunk.loaded) continue;
            //if (chunk.frame >= this.frame - this.ageLimit) continue; // not old
            // - explicit position check rather than age:
            const inBounds = p => {
                for (let k of ['x', 'y', 'z']) {
                    if (p[k] < minChunkPos[k] ||
                        p[k] > maxChunkPos[k]) {
                        return false;
                    }
                }
                return true;
            }
            if (inBounds(chunk.id)) continue;
            
            delete this.posToChunk[chunkPos];
            if (chunk.mesh) {
                this.scene.remove(chunk.mesh);
                chunk.mesh.geometry.dispose();
            }
            ++this.numRemoved;
        }
        
        {
            const elapsed = Date.now() - t1;
            //if (elapsed) console.log('Update in ' + elapsed);
        }
        
        this.culledDistant = numDistant;
        this.numVisited = numVisited;
        if (false) // - debugging
        {
            if (!this.dontShowTooOften) this.dontShowTooOften = 1;
            if (this.dontShowTooOften % 60 == 0)
            console.log(`${numVisited}/${this.chunkOrder.length} chunks updated in ${t2 - t1}, removal in ${Date.now() - t2}`);
            ++this.dontShowTooOften;
        }
        
        ++this.frame;
    }
    
    // f: function with parameters (bi, {x: ix, y: iy, z: iz})
    blockIteration(chunk, f) {
        let bi = 0;
        let c = new THREE.Vector3();
        for (c.z = 0; c.z < this.chunkResolution.z; ++c.z)
        for (c.y = 0; c.y < this.chunkResolution.y; ++c.y)
        for (c.x = 0; c.x < this.chunkResolution.x; ++c.x) {
            f(bi, c);
            ++bi;
        }
    }
}
