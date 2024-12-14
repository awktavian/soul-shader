// Quantum Time-Crystal Shader (Vermeer-Refined)
//
// This GLSL code can be directly pasted into Shadertoy. To use in other contexts:
// - Provide uniforms: iTime (float), iResolution (vec3).
// - Implement a camera and ray direction logic as in mainImage.
// - Integrate the raymarching into your rendering pipeline.
//
// For Shadertoy: Ensure iTime and iResolution are set. Remove or ignore any sections
// not needed in your engine. The code is self-contained for Shadertoy.
//
// Artistic parameters and their emotional impact explained in README.md.

#ifdef GL_ES
precision mediump float;
#endif

uniform float iTime;
uniform vec3  iResolution;

#define MAX_STEPS 128
#define MAX_DIST 20.0
#define SURF_DIST 0.0005
#define VOL_STEPS 20
#define VOL_STEP_SIZE 0.05
#define CUBE_SIZE 10.0
#define MAX_BOUNCES 5

float fRefrIndex          = 1.38;   
float fBaseDispersion     = 0.015;  
float fAtomRadius         = 0.16;   
float fLatticePeriod      = 4.5;    
float fNoiseScale         = 2.5;    
float fVolScatStrength    = 0.5;    
float fVolAbsorption      = 0.18;   
float fMorphSpeed         = 0.4;    
float fStarGlowIntensity  = 40.0;

float gGlitch = 0.0;

float hash21(vec2 p) {
    p = fract(p*vec2(123.34,234.45));
    p += dot(p,p+45.32);
    return fract(p.x*p.y);
}

float noise3D(vec3 p) {
    vec3 i = floor(p);
    vec3 f = fract(p);
    f = f*f*(3.0-2.0*f);

    float n000 = hash21(i.xy);
    float n100 = hash21((i+vec3(1,0,0)).xy);
    float n010 = hash21((i+vec3(0,1,0)).xy);
    float n110 = hash21((i+vec3(1,1,0)).xy);
    float n001 = hash21((i+vec3(0,0,1)).xy);
    float n101 = hash21((i+vec3(1,0,1)).xy);
    float n011 = hash21((i+vec3(0,1,1)).xy);
    float n111 = hash21((i+vec3(1,1,1)).xy);

    float nx00 = mix(n000,n100,f.x);
    float nx10 = mix(n010,n110,f.x);
    float nx01 = mix(n001,n101,f.x);
    float nx11 = mix(n011,n111,f.x);
    float nxy0 = mix(nx00,nx10,f.y);
    float nxy1 = mix(nx01,nx11,f.y);
    float nxyz = mix(nxy0,nxy1,f.z);
    return nxyz;
}

mat3 rotateY(float a) {
    float s=sin(a), c=cos(a);
    return mat3(c,0,-s, 0,1,0, s,0,c);
}

vec3 hsv2rgb(vec3 c) {
    vec3 K=vec3(1.0,2.0/3.0,1.0/3.0);
    vec3 p=abs(fract(c.xxx+K.xyz)*6.0-3.0);
    return c.z*mix(vec3(1.0),clamp(p-1.0,0.0,1.0),c.y);
}

float fbm(vec3 p) {
    float f=0.0;
    float amp=0.5;
    for(int i=0; i<5; i++) {
        f += amp*noise3D(p);
        p *= 2.01;
        amp *= 0.5;
    }
    return f;
}

vec3 adjustForGlitchVermeer(vec3 hsvColor) {
    if(gGlitch > 0.0) {
        float lavenderHue = 0.7;
        hsvColor.x = mix(hsvColor.x, lavenderHue, 0.2*gGlitch); 
        hsvColor.y = mix(hsvColor.y, 0.8*hsvColor.y, gGlitch);
        hsvColor.z += 0.03*gGlitch;
    }
    return hsvColor;
}

vec3 auroraBackground(vec3 rd, float t) {
    vec3 p = vec3(rd.x*2.0, rd.z*2.0, t*0.1); 
    float base = fbm(p*1.5);
    float streaks = fbm(p + vec3(0.0,10.0,0.0)*0.5);
    float pattern = smoothstep(0.3,1.0, base + streaks*0.5);

    float hue = 0.4 + 0.08*sin(t*0.2 + rd.z*2.0);
    float sat = 0.9;
    float val = 0.4 + 0.4*pattern;
    vec3 hsvAurora = vec3(hue,sat,val);
    hsvAurora = adjustForGlitchVermeer(hsvAurora);
    vec3 auroraColor = hsv2rgb(hsvAurora);

    float foam = fbm(p*4.0 + vec3(t*0.5, t*0.3, t*0.7));
    foam = smoothstep(0.4,1.0, foam)*0.25;
    vec3 foamHsv = vec3(hue+0.2,0.5,0.6);
    foamHsv = adjustForGlitchVermeer(foamHsv);
    vec3 foamColor = hsv2rgb(foamHsv)*foam;

    vec3 color = auroraColor + foamColor;
    float horizonFade = smoothstep(-0.3,0.0,rd.y);
    color *= horizonFade;
    return color;
}

const vec3 caffeineSubs[3] = vec3[3](
    vec3(1.2,0.3,0.0),
    vec3(-0.8,0.4,1.0),
    vec3(0.0,0.3,-1.2)
);

const vec3 theobromineSubs[2] = vec3[2](
    vec3(1.2,0.3,0.0),
    vec3(-0.8,0.4,1.0)
);

const float ringOffsetCaffeine[8] = float[8](0.0,0.1,-0.1,0.0,0.05,-0.05,0.1,-0.1);
const float ringOffsetTheo[8]     = float[8](0.1,0.0,0.0,0.1,-0.05,0.05,0.0,0.0);

vec3 ringAtom(float i, float N) {
    float angle = 2.0*3.14159*i/N;
    return vec3(0.7*cos(angle),0.0,0.7*sin(angle));
}

void getMoleculePositions(out vec3 atoms[16], float morph) {
    for (int i=0; i<8; i++) {
        float yC = ringOffsetCaffeine[i];
        float yT = ringOffsetTheo[i];
        float y = mix(yC, yT, morph);
        vec3 basePos = ringAtom(float(i),8.0) + vec3(0.0,y,0.0);
        atoms[i] = basePos;
    }

    for (int i=0; i<3; i++) {
        vec3 cpos = caffeineSubs[i];
        vec3 tpos = (i<2) ? theobromineSubs[i] : vec3(0.0);
        float w = (i<2) ? 1.0 : (1.0 - morph);
        vec3 pos = mix(cpos, tpos, morph)*w;
        atoms[8+i] = pos;
    }

    for (int i=11; i<16; i++) {
        atoms[i] = vec3(0.0);
    }
}

float moleculeSDF(vec3 p, float morph) {
    vec3 atoms[16];
    getMoleculePositions(atoms, morph);
    float scale = 1.5 + 0.1*sin(morph*3.14159);
    mat3 rot = rotateY(morph*3.14159);
    p = rot * p * scale;

    float d = 999.0;
    for (int i=0; i<16; i++) {
        if (i>=11 && atoms[i]==vec3(0.0)) break;
        float distToAtom = length(p - rot*atoms[i]) - fAtomRadius;
        d = min(d, distToAtom);
    }
    return d;
}

float quantumMorph(float t) {
    float baseMorph = 0.5 + 0.5*sin(t*fMorphSpeed*2.0);
    float glitchFreq = 0.3;
    float g = fract(t*glitchFreq);
    float glitch = smoothstep(0.0,0.03,abs(g-0.1)) * 0.6;
    gGlitch = glitch;
    baseMorph = mix(baseMorph, 1.0-baseMorph, glitch);
    baseMorph += (noise3D(vec3(t*0.5, t*0.3, t*0.7)) - 0.5)*0.1;
    return clamp(baseMorph,0.0,1.0);
}

float timeCrystalSDF(vec3 p, float t) {
    float morph = quantumMorph(t);
    vec3 cell = mod(p, fLatticePeriod) - 0.5*fLatticePeriod;
    float baseDist = moleculeSDF(cell, morph);
    float n = (noise3D(p*fNoiseScale + vec3(t*0.3, t*0.2, t*0.4))-0.5)*0.1;
    baseDist += n;
    return baseDist;
}

float raymarchSurface(vec3 ro, vec3 rd, float t) {
    float dist = 0.0;
    for (int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd*dist;
        float d = timeCrystalSDF(p,t);
        if (d < SURF_DIST) return dist;
        dist += d;
        if (dist > MAX_DIST) break;
    }
    return -1.0;
}

vec3 getNormal(vec3 p, float t) {
    float d = timeCrystalSDF(p,t);
    vec2 e = vec2(0.001,0);
    vec3 n = vec3(
        timeCrystalSDF(p+vec3(e.x,e.y,e.y), t)-d,
        timeCrystalSDF(p+vec3(e.y,e.x,e.y), t)-d,
        timeCrystalSDF(p+vec3(e.y,e.y,e.x), t)-d
    );
    return normalize(n);
}

float getDispersion(float t) {
    return fBaseDispersion * (1.0 + 0.3*sin(t*0.7));
}

float fresnel(vec3 n, vec3 v, float ior) {
    float f0 = pow((ior-1.0)/(ior+1.0),2.0);
    float cosi = dot(n,-v);
    return f0 + (1.0 - f0)*pow(1.0 - cosi,5.0);
}

vec3 refrDirR(vec3 I, vec3 N, float disp){return refract(I,N,1.0/(fRefrIndex - disp));}
vec3 refrDirG(vec3 I, vec3 N, float disp){return refract(I,N,1.0/fRefrIndex);}
vec3 refrDirB(vec3 I, vec3 N, float disp){return refract(I,N,1.0/(fRefrIndex + disp));}

float raymarchExit(vec3 ro, vec3 rd, float t) {
    float dist = 0.0;
    for (int i=0; i<MAX_STEPS; i++) {
        vec3 p = ro + rd*dist;
        float d = timeCrystalSDF(p,t);
        if (d > 0.1) return dist; 
        dist += 0.05;
        if (dist>MAX_DIST) break;
    }
    return -1.0;
}

vec3 quantumScatteringColor(vec3 pos, float t) {
    float hue = fract(pos.x*0.2 + pos.y*0.3 + pos.z*0.2 + t*0.5);
    float sat = 0.8 + 0.2*sin(t+pos.y);
    float val = 0.45+0.5*sin(pos.x*5.0+pos.z*5.0+t*1.0);
    vec3 hsv = vec3(hue, clamp(sat,0.0,1.0), clamp(val,0.0,1.0));
    hsv = adjustForGlitchVermeer(hsv);
    return hsv2rgb(hsv);
}

bool intersectCube(vec3 ro, vec3 rd, float size, out float distHit, out vec3 normal) {
    vec3 invDir = 1.0/rd;
    vec3 tminv = (-vec3(size)-ro)*invDir;
    vec3 tmaxv = ( vec3(size)-ro)*invDir;

    vec3 tmin = min(tminv,tmaxv);
    vec3 tmax = max(tminv,tmaxv);

    float tNear = max(max(tmin.x,tmin.y),tmin.z);
    float tFar  = min(min(tmax.x,tmax.y),tmax.z);

    if(tFar<0.0 || tNear>tFar) return false;

    distHit = tNear>0.0 ? tNear : tFar;

    vec3 hitPos = ro+rd*distHit;
    normal = vec3(0.0);
    float eps = 0.001;
    if(abs(hitPos.x + size)<eps) normal=vec3(-1,0,0);
    else if(abs(hitPos.x - size)<eps) normal=vec3(1,0,0);
    else if(abs(hitPos.y - size)<eps) normal=vec3(0,1,0);
    else if(abs(hitPos.y + size)<eps) normal=vec3(0,-1,0);
    else if(abs(hitPos.z - size)<eps) normal=vec3(0,0,1);
    else if(abs(hitPos.z + size)<eps) normal=vec3(0,0,-1);

    return true;
}

vec3 cubeOrbsColor(vec3 p, float t) {
    vec3 orbPos[3] = vec3[3](
        vec3(2.0,1.0,0.0),
        vec3(-1.5,-1.0,2.0),
        vec3(0.0,1.5,-2.0)
    );
    vec3 col = vec3(0.0);
    for(int i=0;i<3;i++){
        float d = length(p-orbPos[i]);
        float intensity = exp(-4.0*d);
        float hue = fract(t*0.2 + float(i)*0.3);
        vec3 hsv = vec3(hue,0.9,1.0);
        hsv = adjustForGlitchVermeer(hsv);
        col += hsv2rgb(hsv)*intensity;
    }
    return col;
}

vec3 integrateVolume(vec3 startPos, vec3 rd, float t) {
    vec3 pos = startPos + rd*0.01; 
    float exitDist = raymarchExit(pos, rd, t);
    if (exitDist<0.0) exitDist = 2.0;
    float stepSize = exitDist/float(VOL_STEPS);
    vec3 col = vec3(0.0);
    float transmittance = 1.0;

    for (int i=0; i<VOL_STEPS; i++) {
        vec3 samplePos = pos + rd*(float(i)+0.5)*stepSize;
        float density = 0.5 + 0.5*noise3D(samplePos*fNoiseScale + vec3(t));
        density = smoothstep(0.3,1.0,density)*0.6;

        vec3 scatColor = quantumScatteringColor(samplePos, t);
        float dMol = timeCrystalSDF(samplePos, t);
        float absorb = exp(-fVolAbsorption*density*stepSize);
        float scatter = density * fVolScatStrength*stepSize;

        vec3 starGlow = vec3(0.0);
        if (dMol < 0.03) {
            float coolFactor = 1.0 + 0.3*gGlitch; 
            starGlow = vec3(1.1*coolFactor,1.0,0.85*coolFactor)*fStarGlowIntensity*exp(-20.0*dMol);
        }

        col += transmittance * (scatter * scatColor + starGlow * stepSize);
        transmittance *= absorb;
        if (transmittance<0.001) break;
    }

    col += auroraBackground(rd,t)*transmittance*0.2;
    return col;
}

// mainImage for Shadertoy. For other apps, adapt camera & uniforms.
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = (fragCoord - 0.5*iResolution.xy)/iResolution.y;
    float t = iTime;

    vec3 ro = vec3(0.0, 0.0, 5.0);
    vec3 ta = vec3(0.0,0.0,0.0);
    vec3 ww = normalize(ta - ro);
    vec3 uu = normalize(cross(vec3(0,1,0), ww));
    vec3 vv = cross(ww, uu);
    vec3 rd = normalize(uv.x*uu + uv.y*vv + ww);
    rd = rotateY(t*0.2)*rd;

    float d = raymarchSurface(ro, rd, t);
    vec3 col;
    if (d>0.0) {
        vec3 p = ro + rd*d;
        vec3 n = getNormal(p, t);
        float refl = fresnel(n, rd, fRefrIndex);

        float disp = getDispersion(t);
        vec3 pInside = p - n*0.001;
        vec3 volR = integrateVolume(pInside, refrDirR(rd,n,disp), t);
        vec3 volG = integrateVolume(pInside, refrDirG(rd,n,disp), t);
        vec3 volB = integrateVolume(pInside, refrDirB(rd,n,disp), t);
        vec3 refrCol = vec3(volR.r, volG.g, volB.b);

        vec3 reflDir = reflect(rd,n);
        vec3 reflOrigin = p + n*0.001;
        vec3 reflCol = infinityMirrorEnv(reflOrigin, reflDir, t);

        col = mix(refrCol, reflCol, refl);
    } else {
        col = auroraBackground(rd, t)*0.8;
    }

    col = col/(col+vec3(1.0));
    col = pow(col, vec3(0.4545));
    fragColor = vec4(col,1.0);
}
