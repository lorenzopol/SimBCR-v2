#version 330
struct Light {
    vec3 position;
    vec3 Ia;
    vec3 Id;
    vec3 Is;
};

layout (location = 0) out vec4 fragColor;

in vec2 uv_0;
in vec3 normal;
in vec3 fragPos;


uniform sampler2D u_texture_0;
uniform Light light;
uniform vec3 camPos;


vec3 getLight(vec3 color) {
    vec3 Normal = normalize(normal);
    vec3 ambient = light.Ia;

    // diffuse light (Lambert law)
    vec3 lightDir = normalize(light.position - fragPos);
    float diff  = max(0, dot(lightDir, Normal));
    vec3 diffuse = diff * light.Id;

    // specular
    vec3 viewDir = normalize(camPos - fragPos);
    vec3 reflectDir = reflect(-lightDir, Normal);
    // values bigger than 32 will make the surface shinier, less then 32 will make it roughier
    float spec = pow(max(dot(viewDir, reflectDir), 0), 32);
    vec3 specular = spec * light.Is;

    return color * (ambient + diffuse + specular);
}
void main() {
    float gamma = 2.2;

    vec3 color = texture(u_texture_0, uv_0).rgb;
    // 1. undo gamma correction (done in the tex)
    color = pow(color, vec3(gamma));

    // 2. calculate light
    color = getLight(color);

    // 3. gamma correct the blended light + tex
    color = pow(color, 1/vec3(gamma));
    fragColor = vec4(color, 1.0);
}