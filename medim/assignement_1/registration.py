import cv2 as cv2
import numpy as np
from os import listdir

def tranform_rigid_from_matches(matches):
    x, y = zip(*matches)

    x_bar = np.sum(x, axis=0) / len(x)
    y_bar = np.sum(y, axis=0) / len(y)

    x_tilde = x - x_bar
    y_tilde = y - y_bar

    H = np.einsum('bi,bo->bio', x_tilde, y_tilde)
    H = np.sum(H, axis=0)
    u, _, vh = np.linalg.svd(H, full_matrices=False)

    r_star = np.matmul(u, vh)
    if (np.linalg.det(r_star) < 0):
        r_star[:1] = -r_star[:1]
    translation = y_bar - np.matmul(r_star, x_bar)
    return r_star, translation


def sqaured_distance_of_matches(keypoint_pairs):
    sq_dist = 0
    for kp1, kp2 in keypoint_pairs:
        sq_dist += np.sum(np.square(kp1 - kp2))
    return sq_dist / len(keypoint_pairs)

def draw_keypoints(image, keypoints):
    for kp in keypoints:
        x = int(np.round(kp[0]))
        y = int(np.round(kp[1]))
        try:
            image[x][y] = 0
            image[x - 1][y + 1] = 0
            image[x + 1][y + 1] = 0
            image[x - 1][y - 1] = 0
            image[x + 1][y - 1] = 0
        except IndexError:
            pass
    return image
def draw_keypoints2(image, keypoints):
    for kp in keypoints:
        x = int(np.round(kp[0]))
        y = int(np.round(kp[1]))
        try:
            image[x, y] = 0
            image[x, y + 1] = 0
            image[x, y - 1] = 0
            image[x - 1, y] = 0
            image[x + 1, y] = 0
        except IndexError:
            pass
    return image

def get_top_matches(kp1, kp2, matches, number_of_matches_to_keep = 50):
    sorted_matches = sorted(matches, key=lambda x: x.distance)
    top_matches = sorted_matches[0:1] = sorted_matches[2:5]
    # top_matches = sorted_matches[0:5]
    top_kp1 = [kp1[m.queryIdx] for m in top_matches]
    top_kp2 = [kp2[m.trainIdx] for m in top_matches]
    return top_kp1, top_kp2, top_matches

def show_images(size, images, axis=1):
    vis = np.concatenate(
        tuple(cv2.resize(img, size) for img in images), axis=axis)
    cv2.imshow("boop", vis)
    cv2.waitKey(0)
    cv2.destroyAllWindows()

def is_within_distance(x, y, distance):
    return np.sum(np.square(x - y)) < distance

def get_filenames():
    return [file_name for file_name in listdir('./data/coll_1/HE/') if file_name[-3:] == 'bmp']
    # return ['8.0.bmp', '8.3.bmp', '8.5.bmp', '1.3.bmp', '1.5.bmp', '6.6.bmp', '6.7.bmp']#[5:6]
    # return ['8.4.bmp']

def get_image_pair(file_name):
    return cv2.imread("./data/coll_1/HE/" + file_name), cv2.imread("./data/coll_1/p63AMACR/" + file_name)

def down_sample_image(img, scale = 0.5):
    return cv2.resize(img, (0, 0), fx=scale, fy=scale)

def greyscale_image(img):
    return cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

def get_sift_matches(img1, img2):
    sift = cv2.xfeatures2d.SIFT_create()

    kp1, des1 = sift.detectAndCompute(img1, None)
    kp2, des2 = sift.detectAndCompute(img2, None)

    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
    matches = bf.match(des1, des2)
    nbr_of_top_matches = 20# int(len(matches) / 10)
    matches = sorted(matches, key = lambda x:x.distance)[:nbr_of_top_matches]

    # return [tuple([kp1.pt, kp2.pt]) for kp1, kp2 in zip(kp1, kp2)]
    return [tuple([kp1[m.queryIdx].pt, kp2[m.trainIdx].pt]) for m in matches]

def random_sample(list, sample_size):
    indices = np.random.choice(len(list)), 2
    return [list[i] for i in indices]

def rotate_points(points, R, t):
    return np.array([np.matmul(R, p) + t for p in points])

def points_within_treshhold(x, y, treshhold):
    return [i for i, p in enumerate(np.sum(np.square(np.array(y)-np.array(x)), axis=1)) if p < treshhold]

def transform_image(image, R, t):
    rows, cols, _ = image.shape
    center_x = rows / 2
    center_y = cols / 2
    M = np.array([[R[0, 0], R[0, 1], (1 - R[0, 0]) * (center_x + t[0]) - R[0, 1] * (center_y+t[1])],
                   [R[1, 0], R[1, 1], R[0, 1] * (center_x + t[0]) + (1 - R[0, 0]) * (center_y+t[1])]])
    return cv2.warpAffine(image, M.astype(np.float32), (cols, rows))

def transform_image3(image, R, t):
    rows, cols, _ = image.shape
    center_x = cols / 2
    center_y = rows / 2
    M = np.array([[R[0, 0], R[0, 1], (1 - R[0, 0]) * center_x - R[0, 1] * center_y],
                  [R[1, 0], R[1, 1], R[0, 1] * center_x + (1 - R[0, 0]) * center_y]])
    print('RM')
    print(R)
    print(M)
    M2 = np.array([[1, 0, t[0]],
                   [0, 1, t[1]]])
    rotated = cv2.warpAffine(image, M.astype(np.float32), (cols, rows))
    return cv2.warpAffine(rotated, M2.astype(np.float32), (cols, rows))


def transform_image2(image, R, t):
    rows, cols, _ = image.shape
    center_x = rows / 2
    center_y = cols / 2
    M = np.array([[R[0, 0], R[0, 1], (1 - R[0, 0]) * center_x - R[0, 1] * center_y],
                   [R[1, 0], R[1, 1], R[0, 1] * center_x + (1 - R[0, 0]) * center_y]])
    rows, cols, _ = image.shape
    return cv2.warpAffine(image, M.astype(np.float32), (cols, rows))

def overlay_images(size, img1, img2):
    i1 = cv2.resize(img1, size)
    i2 = cv2.resize(img2, size)
    i1[:, ::2] = 0
    i2[:, 1::2] = 0
    return i1 + i2

if __name__ == "__main__":
    for file_name in get_filenames():
        img1, img2 = get_image_pair(file_name)
        scaling = 0.4
        small1 = down_sample_image(img1, scaling)
        small2 = down_sample_image(img2, scaling)
        gray1 = greyscale_image(small1)
        gray2 = greyscale_image(small2)
        matches = get_sift_matches(gray1, gray2)
        best_match = 0
        best_points = None
        for i in range(300):
            matches_sample = random_sample(matches, 3)
            r_star, translation = tranform_rigid_from_matches(matches_sample)
            x_rotated = rotate_points([x for x, _ in matches], r_star, translation)
            y = np.array([y for _, y in matches])
            treshhold = 50 #gray1.shape[0] / 10
            good_point_indices = points_within_treshhold(x_rotated, y, treshhold)
            if len(good_point_indices) > best_match:
                best_match = len(good_point_indices)
                best_points = good_point_indices
        good_points = [matches[i] for i in best_points]
        print(file_name)
        print(len(good_points))
        r_star, translation = tranform_rigid_from_matches(good_points)
        print('t')
        print(translation)
        rotated_img1 = transform_image3(small1, r_star, translation)

        overlayed = overlay_images((400, 400), rotated_img1, small2)

        rows, cols, _ = small1.shape
        rotated_keypoints = [np.matmul(r_star, kp1) + translation for kp1, _ in good_points]
        kp2 = [kp2 for _, kp2 in good_points]
        kp1 = [kp1 for kp1, _ in good_points]
        white_kp1_rotated = draw_keypoints2(np.full((rows, cols, 3), 255, np.uint8), rotated_keypoints)
        white_both = draw_keypoints(white_kp1_rotated, kp2)

        white_kp1_clean = draw_keypoints2(np.full((rows, cols, 3), 255, np.uint8), kp1)
        white_both_clean = draw_keypoints(white_kp1_clean, kp2)

        show_images((400, 400), (white_both_clean, white_both, overlayed))

if __name__ != "__main__":
    image1 = cv2.imread("./data/coll_1/HE/1.1.bmp")
    image2 = cv2.imread("./data/coll_1/p63AMACR/1.1.bmp")

    small1 = cv2.resize(image1, (0, 0), fx=0.2, fy=0.2)
    small2 = cv2.resize(image2, (0, 0), fx=0.2, fy=0.2)

    gray1 = cv2.cvtColor(small1, cv2.COLOR_BGR2GRAY)
    gray2 = cv2.cvtColor(small2, cv2.COLOR_BGR2GRAY)

    sift = cv2.xfeatures2d.SIFT_create()

    kp1, des1 = sift.detectAndCompute(gray1, None)
    kp2, des2 = sift.detectAndCompute(gray2, None)

    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
    matches = bf.match(des1, des2)

    top_kp1, top_kp2, top_matches = get_top_matches(kp1, kp2, matches, 5)

    show_images((600, 600), (draw_keypoints(small1, top_kp1), draw_keypoints(small2, top_kp2)))

    r_star, translation = tranform_rigid_from_matches(top_matches)

    num_rows, num_cols, _ = small1.shape
    translation_x = (1 - r_star[0][0]) * num_cols / 2 - (r_star[0][1]) * num_rows / 2
    translation_y = (r_star[0][1]) * num_cols / 2 + (1 - r_star[0][0]) * num_rows / 2


    M = np.hstack((r_star, np.array([[translation[0]], [translation[1]]])))

    downscale_size = (600, 600)
    rotated = cv2.warpAffine(small1, M.astype(np.float32), (num_cols, num_rows))
    rotated_keypoints = [np.matmul(r_star, kp.pt) + translation for kp in top_kp1]
    rotated_withkeypoints = draw_keypoints2(rotated, rotated_keypoints)
    show_images((600, 600), (rotated_withkeypoints, draw_keypoints(small2, top_kp2)))

    white_kp1_rotated = draw_keypoints2(np.full((num_rows, num_cols, 3), 255, np.uint8), rotated_keypoints)
    white_both = draw_keypoints(white_kp1_rotated, top_kp2)
    white_kp2 = draw_keypoints(np.full((num_rows, num_cols, 3), 255, np.uint8), top_kp2)

    white_kp1 = draw_keypoints2(np.full((num_rows, num_cols, 3), 255, np.uint8), [kp.pt for kp in top_kp1])
    both_no_rotate = draw_keypoints(white_kp1, top_kp2)
    show_images((600, 600), [both_no_rotate, white_both])

    kp1pts = [np.array(kp1.pt) for kp1 in top_kp1]
    kp2pts = [np.array(kp2.pt) for kp2 in top_kp2]
    rotkpt = [rotated_kp for rotated_kp in rotated_keypoints]

    keypoint_pairs = [(kp1, kp2) for kp1, kp2 in zip(kp1pts, kp2pts)]
    rotated_keypoint_pairs = [(kp1, kp2) for kp1, kp2 in zip(rotkpt, kp2pts)]

    sq = sqaured_distance_of_matches(keypoint_pairs)
    print(sq)
    sq = sqaured_distance_of_matches(rotated_keypoint_pairs)
    print(sq)
